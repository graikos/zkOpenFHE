#ifndef OPENFHE_PROOFSYSTEM_LIBSNARK_CPP
#define OPENFHE_PROOFSYSTEM_LIBSNARK_CPP

#include "proofsystem/proofsystem_libsnark.h"
#include "proofsystem/gadgets_libsnark.h"
#include "schemerns/rns-cryptoparameters.h"
#include "keyswitch/keyswitch-bv.h"
#include "scheme/bgvrns/bgvrns-leveledshe.h"

#include <vector>
#undef NDEBUG
#include <cassert>

#define LIBSNARK_PROOF_METADATA_KEY "libsnark_proof_metadata"

using namespace libsnark;
using std::cout, std::endl;
using std::vector;

FieldT get_max_field_element(const vector<FieldT>& vec) {
    auto it = std::max_element(vec.begin(), vec.end(), [](const FieldT& x, const FieldT& y) { return lt(x, y); });
    return *it;
}

std::shared_ptr<LibsnarkProofMetadata> LibsnarkProofSystem::GetProofMetadata(const Ciphertext<DCRTPoly>& ciphertext) {
    auto it = ciphertext->FindMetadataByKey(LIBSNARK_PROOF_METADATA_KEY);
    if (ciphertext->MetadataFound(it)) {
        return std::dynamic_pointer_cast<LibsnarkProofMetadata>(ciphertext->GetMetadata(it));
    }
    else {
        OPENFHE_THROW(openfhe_error, "Attempt to access metadata (LibsnarkProofMetadata) that has not been set.");
    }
}

void LibsnarkProofSystem::SetProofMetadata(const Ciphertext<DCRTPoly>& ciphertext,
                                           const std::shared_ptr<LibsnarkProofMetadata>& metadata) {
    ciphertext->SetMetadataByKey(LIBSNARK_PROOF_METADATA_KEY, metadata);
}

void LibsnarkProofSystem::ConstrainPublicInput(Ciphertext<DCRTPoly>& ciphertext) {
    const size_t num_polys = ciphertext->GetElements().size();
    const size_t num_limbs = ciphertext->GetElements()[0].GetNumOfElements();

    LibsnarkProofMetadata out(num_polys);
    out.max_value = vector<vector<FieldT>>(num_polys);
    out.modulus   = vector<size_t>(num_limbs);

    for (size_t j = 0; j < num_limbs; j++) {
        out.modulus[j] = ciphertext->GetElements()[0].GetElementAtIndex(j).GetModulus().ConvertToInt<unsigned long>();
    }

    for (size_t i = 0; i < ciphertext->GetElements().size(); i++) {
        const auto c_i   = ciphertext->GetElements()[i];
        out[i]           = vector<vector<pb_linear_combination<FieldT>>>(c_i.GetNumOfElements());
        out.max_value[i] = vector<FieldT>(c_i.GetNumOfElements());
        for (size_t j = 0; j < c_i.GetNumOfElements(); j++) {
            const auto c_ij     = c_i.GetElementAtIndex(j);
            const auto& v_ij    = c_ij.GetValues();
            out[i][j]           = vector<pb_linear_combination<FieldT>>(v_ij.GetLength());
            out.max_value[i][j] = FieldT(c_ij.GetModulus().ConvertToInt<unsigned long>()) - 1;

            for (size_t k = 0; k < v_ij.GetLength(); k++) {
                pb_variable<FieldT> tmp;
                tmp.allocate(pb, ciphertext->SerializedObjectName() + "[" + std::to_string(i) + "][" +
                                     std::to_string(j) + "][" + std::to_string(k) + "] (input)");
                pb.val(tmp)  = v_ij[k].ConvertToInt<unsigned long>();
                out[i][j][k] = pb_linear_combination<FieldT>(tmp);
            }
        }
    }

    pb.set_input_sizes(pb.num_inputs() + (out.size() * out[0].size() * out[0][0].size()));
    SetProofMetadata(ciphertext, std::make_shared<LibsnarkProofMetadata>(out));
}

void LibsnarkProofSystem::constrain_addmod_lazy(const LibsnarkProofMetadata& in1, const size_t index_1,
                                                const LibsnarkProofMetadata& in2, const size_t index_2,
                                                LibsnarkProofMetadata& out, const size_t index_out) {
    assert(index_1 < in1.size());
    assert(index_2 < in2.size());
    assert(index_out < out.size());
    auto num_limbs = in1[index_1].size();
    assert(in2[index_2].size() == num_limbs);
    assert(out[index_out].size() == num_limbs);
    auto modulus = in1.modulus;
    assert(in2.modulus == modulus);
    assert(out.modulus == modulus && "modulus of `out' is not set");

    vector<FieldT> out_max_value(num_limbs);
    vector<bool> field_overflow(num_limbs);
    out.max_value[index_out] = vector<FieldT>(out[index_out].size());
    for (size_t j = 0; j < num_limbs; ++j) {
        size_t out_bit_size = std::max(in1.get_bit_size(index_1, j), in2.get_bit_size(index_2, j)) + 1ul;
        out_max_value[j]    = in1.max_value[index_1][j] + in2.max_value[index_2][j];
        field_overflow[j]   = out_bit_size >= FieldT::num_bits;

        if (field_overflow[j]) {
            // Eager witness generation, add modulus constraints
            auto g = BatchGadget<FieldT, AddModGadget<FieldT>>(pb, in1[index_1][j], in1.max_value[index_1][j],
                                                               in2[index_2][j], in2.max_value[index_2][j], modulus[j]);
            g.generate_r1cs_constraints();
            g.generate_r1cs_witness();
            out[index_out][j]           = g.get_output();
            out.max_value[index_out][j] = FieldT(modulus[j]) - 1;
        }
        else {
            // Lazy branch, do not add modulus constraints, but track size of values for later
            out[index_out][j] = vector<pb_linear_combination<FieldT>>(in1[index_1][j].size());
            for (size_t k = 0; k < out[index_out][j].size(); ++k) {
                out[index_out][j][k].assign(pb, in1[index_1][j][k] + in2[index_2][j][k]);
            }
            out.max_value[index_out][j] = out_max_value[j];
        }
    }
}

void LibsnarkProofSystem::constrain_submod_lazy(const LibsnarkProofMetadata& in1, const size_t index_1,
                                                const LibsnarkProofMetadata& in2, const size_t index_2,
                                                LibsnarkProofMetadata& out, const size_t index_out) {
    assert(index_1 < in1.size());
    assert(index_2 < in2.size());
    assert(index_out < out.size());
    auto num_limbs = in1[index_1].size();
    assert(in2[index_2].size() == num_limbs);
    assert(out[index_out].size() == num_limbs);
    const auto modulus(in1.modulus);
    assert(modulus.size() == num_limbs);
    assert(in2.modulus == modulus);
    assert(out.modulus == modulus && "modulus of `out' is not set");

    vector<FieldT> out_max_value(num_limbs);
    vector<bool> field_overflow(num_limbs);
    vector<bool> in2_lt_modulus(num_limbs);

    out.max_value[index_out] = vector<FieldT>(out[index_out].size());
    for (size_t j = 0; j < num_limbs; ++j) {
        assert(modulus.size() == num_limbs);
        const FieldT curr_mod_fieldT = FieldT(modulus[j]);

        out_max_value[j]    = in1.max_value[index_1][j] + curr_mod_fieldT;
        size_t out_bit_size = std::max(in1.get_bit_size(index_1, j), (size_t)ceil(log2(modulus[j]))) + 1ul;
        field_overflow[j]   = out_bit_size >= FieldT::num_bits;
        in2_lt_modulus[j]   = lt(in2.max_value[index_2][j], curr_mod_fieldT);

        auto in2_ij = in2[index_2][j];
        //        auto in2_ij_curr_bit_size = in2.get_bit_size(index_2, j);
        auto in2_ij_max_value = in2.max_value[index_2][j];
        if (!in2_lt_modulus[j]) {
            // We first need to mod-reduce in2[index_2][j][k] before we can compute its negative
            // TODO: is there a way to compute the negative from the lazy/non-reduced value directly?
            auto g_mod =
                BatchGadget<FieldT, ModGadget<FieldT>>(pb, in2[index_2][j], in2.max_value[index_2][j], modulus[j]);
            g_mod.generate_r1cs_constraints();
            g_mod.generate_r1cs_witness();

            in2_ij           = g_mod.get_output();
            in2_ij_max_value = curr_mod_fieldT;
        }

        if (field_overflow[j]) {
            // Eager witness generation, add modulus constraints
            auto g = BatchGadget<FieldT, SubModGadget<FieldT>>(pb, in1[index_1][j], in1.max_value[index_1][j], in2_ij,
                                                               in2_ij_max_value, modulus[j]);
            g.generate_r1cs_constraints();
            g.generate_r1cs_witness();
            out[index_out][j]           = g.get_output();
            out.max_value[index_out][j] = curr_mod_fieldT - 1;
        }
        else {
            // Lazy branch, do not add modulus constraints, but track size of values for later
            out[index_out][j] = vector<pb_linear_combination<FieldT>>(in1[index_1][j].size());
            for (size_t k = 0; k < out[index_out][j].size(); ++k) {
                out[index_out][j][k].assign(pb, in1[index_1][j][k] + curr_mod_fieldT - in2_ij[k]);
            }
            out.max_value[index_out][j] = out_max_value[j];
        }
    }
}

void LibsnarkProofSystem::constrain_mulmod_lazy(const LibsnarkProofMetadata& in1, const size_t index_1,
                                                const LibsnarkProofMetadata& in2, const size_t index_2,
                                                LibsnarkProofMetadata& out, const size_t index_out) {
    assert(index_1 < in1.size());
    assert(index_1 < in1.max_value.size());
    assert(index_2 < in2.size());
    assert(index_2 < in2.max_value.size());
    assert(index_out < out.size());
    assert(index_out < out.max_value.size());
    auto num_limbs = in1[index_1].size();
    assert(in2[index_2].size() == num_limbs);
    auto modulus = in1.modulus;
    assert(in2.modulus == modulus);
    assert(out.modulus == modulus && "modulus of `out' is not set");

    out[index_out].resize(num_limbs);  // assert(out[index_out].size() == num_limbs);

    vector<FieldT> out_max_value(num_limbs);
    vector<bool> field_overflow(num_limbs);
    out.max_value[index_out] = vector<FieldT>(num_limbs);
    for (size_t j = 0; j < num_limbs; ++j) {
        size_t out_bit_size = in1.get_bit_size(index_1, j) + in2.get_bit_size(index_2, j);
        out_max_value[j]    = in1.max_value[index_1][j] * in2.max_value[index_2][j];
        field_overflow[j]   = out_bit_size >= FieldT::num_bits;

        if (field_overflow[j]) {
            // Eager witness generation, add modulus constraints
            auto g = BatchGadget<FieldT, MulModGadget<FieldT>>(pb, in1[index_1][j], in1.max_value[index_1][j],
                                                               in2[index_2][j], in2.max_value[index_2][j], modulus[j]);
            g.generate_r1cs_constraints();
            g.generate_r1cs_witness();
            out[index_out][j]           = g.get_output();
            out.max_value[index_out][j] = FieldT(modulus[j]) - 1;
        }
        else {
            // Lazy branch, only add quadratic constraint for multiplication without mod-reduction
            auto g = BatchGadget<FieldT, MulGadget<FieldT>>(pb, in1[index_1][j], in2[index_2][j]);
            g.generate_r1cs_constraints();
            g.generate_r1cs_witness();
            //            auto x  = g.get_output();
            //            auto oi = out[index_out].at(j) = x;
            out[index_out][j]           = g.get_output();
            out.max_value[index_out][j] = out_max_value[j];
        }
    }
}

void LibsnarkProofSystem::ConstrainAddition(const Ciphertext<DCRTPoly>& ctxt1, const Ciphertext<DCRTPoly>& ctxt2,
                                            Ciphertext<DCRTPoly>& ctxt_out) {
    const auto in1 = *GetProofMetadata(ctxt1);
    const auto in2 = *GetProofMetadata(ctxt2);

    vector<unsigned long> moduli;
    for (const auto& e : ctxt_out->GetElements()[0].GetAllElements()) {
        moduli.push_back(e.GetModulus().ConvertToInt());
    }
    assert(in1.size() == in2.size());
    assert(in1.modulus == in2.modulus);
    if (in1.modulus != moduli) {
        cout << "in1: " << in1.modulus << endl;
        cout << "in2: " << in2.modulus << endl;
        cout << "moduli: " << moduli << endl;
    }
    assert(in1.modulus == moduli);

    LibsnarkProofMetadata out(in1.size());
    out.max_value = vector<vector<FieldT>>(in1.size());
    out.modulus   = in1.modulus;
    for (size_t i = 0; i < in1.size(); i++) {
        out[i] = vector<vector<pb_linear_combination<FieldT>>>(in1[i].size());
        constrain_addmod_lazy(in1, i, in2, i, out, i);
    }
    SetProofMetadata(ctxt_out, std::make_shared<LibsnarkProofMetadata>(out));
}

void LibsnarkProofSystem::ConstrainAddition(const Ciphertext<DCRTPoly>& ctxt, const Plaintext& ptxt,
                                            Ciphertext<DCRTPoly>& ctxt_out) {
    const auto in1 = *GetProofMetadata(ctxt);
    auto pt        = ptxt->GetElement<DCRTPoly>();
    assert(in1[0].size() == pt.GetNumOfElements());
    assert(in1[0][0].size() == pt.GetLength());

    // CAUTION: ptxt->GetLength() is the number of set slots, not the ring dimension! Use pt.GetLength() instead.
    vector<vector<pb_linear_combination<FieldT>>> in2_lc(in1[0].size());
    vector<FieldT> in2_max_value(in1[0].size());
    const size_t ptxt_modulus = pt.GetElementAtIndex(0).GetModulus().ConvertToInt();
    for (size_t i = 0; i < in2_lc.size(); ++i) {
        in2_lc[i].resize(pt.GetLength());
        for (size_t j = 0; j < pt.GetLength(); ++j) {
            pb_variable<FieldT> var;
            var.allocate(pb, "in2_" + std::to_string(i) + "_" + std::to_string(j));
            // TODO: can we re-use some of the range checks for all entries in the input?
            pb.val(var) = FieldT(pt.GetElementAtIndex(0).GetValues()[j].ConvertToInt());

            // Set max_value to be 1 larger than expected max value to trigger mod reduction
            less_than_constant_gadget<FieldT> g(pb, var, FieldT(ptxt_modulus).as_bigint().num_bits(),
                                                FieldT(ptxt_modulus));
            g.generate_r1cs_constraints();
            g.generate_r1cs_witness();
            in2_lc[i][j] = var;
        }
        in2_max_value[i] = FieldT(ptxt_modulus) - 1;
    }

    // SetFormat(EVALUATION)
    assert(ptxt->GetElement<DCRTPoly>().GetFormat() == EVALUATION);
    vector<vector<pb_linear_combination<FieldT>>> eval_lc = in2_lc;
    vector<FieldT> eval_max_value                         = in2_max_value;

    //    const Plaintext ptxt_eval(ptxt);
    //    ptxt_eval->SetFormat(EVALUATION);
    //    DCRTPoly pt_eval = ptxt_eval->GetElement<DCRTPoly>();
    //    vector<vector<pb_linear_combination<FieldT>>> eval_lc;
    //    vector<FieldT> eval_max_value;
    //    ConstrainSetFormat(EVALUATION, pt, pt_eval, in2_lc, in2_max_value, eval_lc, eval_max_value);

    // Add ptxt from 0-th element of ctxt
    LibsnarkProofMetadata out(in1);
    vector<vector<vector<FieldT>>> out_max_values(in1.size());
    out_max_values[0].resize(in1[0].size());
    for (size_t j = 0; j < in1[0].size(); ++j) {
        out_max_values[0][j].resize(in1[0][j].size());
        for (size_t k = 0; k < in1[0][j].size(); ++k) {
            LazyAddModGadget<FieldT> g(pb, in1[0][j][k], in1.max_value[0][j], eval_lc[j][k], eval_max_value[j],
                                       in1.modulus[j]);
            g.generate_r1cs_constraints();
            g.generate_r1cs_witness();
            out[0][j][k]            = g.out;
            out_max_values[0][j][k] = g.out_max_value;
        }
    }

    for (size_t j = 0; j < out[0].size(); ++j) {
        out.max_value[0][j] = get_max_field_element(out_max_values[0][j]);
    }
    SetProofMetadata(ctxt_out, std::make_shared<LibsnarkProofMetadata>(out));
}

void LibsnarkProofSystem::ConstrainSubstraction(const Ciphertext<DCRTPoly>& ctxt1, const Ciphertext<DCRTPoly>& ctxt2,
                                                Ciphertext<DCRTPoly>& ctxt_out) {
    const auto in1 = *GetProofMetadata(ctxt1);
    const auto in2 = *GetProofMetadata(ctxt2);

    vector<unsigned long> moduli;
    for (const auto& e : ctxt_out->GetElements()[0].GetAllElements()) {
        moduli.push_back(e.GetModulus().ConvertToInt());
    }
    assert(in1.size() == in2.size());
    assert(in1.modulus == in2.modulus);
    assert(in1.modulus == moduli);

    LibsnarkProofMetadata out(in1.size());
    out.max_value = vector<vector<FieldT>>(in1.size());
    out.modulus   = in1.modulus;
    for (size_t i = 0; i < in1.size(); i++) {
        out[i] = vector<vector<pb_linear_combination<FieldT>>>(in1[i].size());
        constrain_submod_lazy(in1, i, in2, i, out, i);
    }
    SetProofMetadata(ctxt_out, std::make_shared<LibsnarkProofMetadata>(out));
}

void LibsnarkProofSystem::ConstrainSubstraction(const Ciphertext<DCRTPoly>& ctxt, const Plaintext& ptxt,
                                                Ciphertext<DCRTPoly>& ctxt_out) {
    const auto in1 = *GetProofMetadata(ctxt);
    auto pt        = ptxt->GetElement<DCRTPoly>();
    assert(in1[0].size() == pt.GetNumOfElements());
    assert(in1[0][0].size() == pt.GetLength());

    // CAUTION: ptxt->GetLength() is the number of set slots, not the ring dimension! Use pt.GetLength() instead.
    vector<vector<pb_linear_combination<FieldT>>> in2_lc(in1[0].size());
    vector<FieldT> in2_max_value(in1[0].size());
    const size_t ptxt_modulus = pt.GetElementAtIndex(0).GetModulus().ConvertToInt();
    for (size_t i = 0; i < in2_lc.size(); ++i) {
        in2_lc[i].resize(pt.GetLength());
        for (size_t j = 0; j < pt.GetLength(); ++j) {
            pb_variable<FieldT> var;
            var.allocate(pb, "in2_" + std::to_string(i) + "_" + std::to_string(j));
            // TODO: can we re-use some of the range checks for all entries in the input?
            pb.val(var) = FieldT(pt.GetElementAtIndex(0).GetValues()[j].ConvertToInt());

            // Set max_value to be 1 larger than expected max value to trigger mod reduction
            less_than_constant_gadget<FieldT> g(pb, var, FieldT(ptxt_modulus).as_bigint().num_bits(),
                                                FieldT(ptxt_modulus));
            g.generate_r1cs_constraints();
            g.generate_r1cs_witness();
            in2_lc[i][j] = var;
        }
        in2_max_value[i] = FieldT(ptxt_modulus) - 1;
    }

    // SetFormat(EVALUATION)
    Plaintext ptxt_eval(ptxt);
    ptxt_eval->SetFormat(EVALUATION);
    DCRTPoly pt_eval = ptxt_eval->GetElement<DCRTPoly>();
    vector<vector<pb_linear_combination<FieldT>>> eval_lc;
    vector<FieldT> eval_max_value;
    ConstrainSetFormat(EVALUATION, pt, pt_eval, in2_lc, in2_max_value, eval_lc, eval_max_value);

    // Subtract ptxt from 0-th element of ctxt
    LibsnarkProofMetadata out(in1);
    vector<vector<vector<FieldT>>> out_max_values(in1.size());
    out_max_values[0].resize(in1[0].size());
    for (size_t j = 0; j < in1[0].size(); ++j) {
        out_max_values[0][j].resize(in1[0][j].size());
        for (size_t k = 0; k < in1[0][j].size(); ++k) {
            LazySubModGadget<FieldT> g(pb, in1[0][j][k], in1.max_value[0][j], eval_lc[j][k], eval_max_value[j],
                                       in1.modulus[j]);
            g.generate_r1cs_constraints();
            g.generate_r1cs_witness();
            out[0][j][k]            = g.out;
            out_max_values[0][j][k] = g.out_max_value;
        }
    }

    for (size_t j = 0; j < out[0].size(); ++j) {
        out.max_value[0][j] = get_max_field_element(out_max_values[0][j]);
    }
    SetProofMetadata(ctxt_out, std::make_shared<LibsnarkProofMetadata>(out));
}

void LibsnarkProofSystem::ConstrainSquare2(const Ciphertext<DCRTPoly>& ctxt, Ciphertext<DCRTPoly>& ctxt_out) {
    const auto in = *GetProofMetadata(ctxt);

    const auto num_limbs = in.modulus.size();
    assert(
        ctxt->GetElements()[0].GetNumOfElements() == ctxt_out->GetElements()[0].GetNumOfElements() &&
        "mismatch between number of limbs between ciphertext input and output. Are you using the FIXEDMANUAL scaling technique?");
    assert(in.size() == 3);

    LibsnarkProofMetadata out(5);
    for (size_t i = 0; i < out.size(); i++) {
        out[i]           = vector<vector<pb_linear_combination<FieldT>>>(num_limbs);
        out.max_value[i] = vector<FieldT>(num_limbs);
    }
    out.modulus = in.modulus;

    LibsnarkProofMetadata tmp(1);
    for (size_t i = 0; i < tmp.size(); i++) {
        tmp[i]           = vector<vector<pb_linear_combination<FieldT>>>(num_limbs);
        tmp.max_value[i] = vector<FieldT>(num_limbs);
    }
    tmp.modulus = in.modulus;

    constrain_mulmod_lazy(in, 0, in, 0, out, 0);

    constrain_mulmod_lazy(in, 0, in, 1, out, 1);
    constrain_addmod_lazy(out, 1, out, 1, out, 1);

    constrain_mulmod_lazy(in, 0, in, 2, tmp, 0);
    constrain_addmod_lazy(tmp, 0, tmp, 0, tmp, 0);
    constrain_mulmod_lazy(in, 1, in, 1, out, 2);
    constrain_addmod_lazy(tmp, 0, out, 2, out, 2);

    constrain_mulmod_lazy(in, 1, in, 2, out, 3);
    constrain_addmod_lazy(out, 3, out, 3, out, 3);

    constrain_mulmod_lazy(in, 2, in, 2, out, 4);

    SetProofMetadata(ctxt_out, std::make_shared<LibsnarkProofMetadata>(out));
}

void LibsnarkProofSystem::ConstrainMultiplication(const Ciphertext<DCRTPoly>& ctxt1, const Ciphertext<DCRTPoly>& ctxt2,
                                                  Ciphertext<DCRTPoly>& ctxt_out) {
    const auto in1 = *GetProofMetadata(ctxt1);
    const auto in2 = *GetProofMetadata(ctxt2);

    vector<unsigned long> moduli;
    for (const auto& e : ctxt_out->GetElements()[0].GetAllElements()) {
        moduli.push_back(e.GetModulus().ConvertToInt());
    }
    assert(in1.size() == in2.size());
    assert(in1.modulus == in2.modulus);
    assert(in1.modulus == moduli);

    assert(in1.size() == 2);
    assert(in2.size() == 2);

    const size_t num_limbs = in1[0].size();

    LibsnarkProofMetadata tmp_0(1);
    LibsnarkProofMetadata tmp_1(1);
    tmp_0[0]           = vector<vector<pb_linear_combination<FieldT>>>(num_limbs);
    tmp_1[0]           = vector<vector<pb_linear_combination<FieldT>>>(num_limbs);
    tmp_0.max_value[0] = vector<FieldT>(num_limbs);
    tmp_1.max_value[0] = vector<FieldT>(num_limbs);
    tmp_0.modulus      = in1.modulus;
    tmp_1.modulus      = in1.modulus;

    LibsnarkProofMetadata out(3);
    out.modulus = in1.modulus;
    for (size_t i = 0; i < out.size(); i++) {
        out[i]           = vector<vector<pb_linear_combination<FieldT>>>(num_limbs);
        out.max_value[i] = vector<FieldT>(num_limbs);
    }
    constrain_mulmod_lazy(in1, 0, in2, 0, out, 0);

    constrain_mulmod_lazy(in1, 0, in2, 1, tmp_0, 0);
    constrain_mulmod_lazy(in1, 1, in2, 0, tmp_1, 0);
    constrain_addmod_lazy(tmp_0, 0, tmp_1, 0, out, 1);

    constrain_mulmod_lazy(in1, 1, in2, 1, out, 2);

    SetProofMetadata(ctxt_out, std::make_shared<LibsnarkProofMetadata>(out));
}

void LibsnarkProofSystem::ConstrainMultiplication(const Ciphertext<DCRTPoly>& ctxt, const Plaintext& ptxt,
                                                  Ciphertext<DCRTPoly>& ctxt_out) {
    const auto in1 = *GetProofMetadata(ctxt);
    auto pt        = ptxt->GetElement<DCRTPoly>();
    assert(in1[0].size() == pt.GetNumOfElements());
    assert(in1[0][0].size() == pt.GetLength());

    // CAUTION: ptxt->GetLength() is the number of set slots, not the ring dimension! Use pt.GetLength() instead.
    vector<vector<pb_linear_combination<FieldT>>> in2_lc(in1[0].size());
    vector<FieldT> in2_max_value(in1[0].size());
    const size_t ptxt_modulus = pt.GetElementAtIndex(0).GetModulus().ConvertToInt();
    for (size_t i = 0; i < in2_lc.size(); ++i) {
        in2_lc[i].resize(pt.GetLength());
        for (size_t j = 0; j < pt.GetLength(); ++j) {
            pb_variable<FieldT> var;
            var.allocate(pb, "in2_" + std::to_string(i) + "_" + std::to_string(j));
            // TODO: can we re-use some of the range checks for all entries in the input?
            pb.val(var) = FieldT(pt.GetElementAtIndex(0).GetValues()[j].ConvertToInt());

            // Set max_value to be 1 larger than expected max value to trigger mod reduction
            less_than_constant_gadget<FieldT> g(pb, var, FieldT(ptxt_modulus).as_bigint().num_bits(),
                                                FieldT(ptxt_modulus));
            g.generate_r1cs_constraints();
            g.generate_r1cs_witness();
            in2_lc[i][j] = var;
        }
        in2_max_value[i] = FieldT(ptxt_modulus) - 1;
    }

    // SetFormat(EVALUATION)
    //    assert(ptxt->GetElement<DCRTPoly>().GetFormat() == EVALUATION);
    const Plaintext ptxt_eval(ptxt);
    ptxt_eval->SetFormat(EVALUATION);
    DCRTPoly pt_eval = ptxt_eval->GetElement<DCRTPoly>();
    vector<vector<pb_linear_combination<FieldT>>> eval_lc;
    vector<FieldT> eval_max_value;
    ConstrainSetFormat(EVALUATION, pt, pt_eval, in2_lc, in2_max_value, eval_lc, eval_max_value);

    // Multiply each entry of ctxt with ptxt
    LibsnarkProofMetadata out(in1.size());
    vector<vector<vector<FieldT>>> out_max_values(in1.size());
    out.modulus = in1.modulus;
    for (size_t i = 0; i < in1.size(); ++i) {
        out[i].resize(in1[i].size());
        out_max_values[i].resize(in1[i].size());
        for (size_t j = 0; j < in1[i].size(); ++j) {
            out[i][j].resize(in1[i][j].size());
            out_max_values[i][j].resize(in1[i][j].size());
            for (size_t k = 0; k < in1[i][j].size(); ++k) {
                LazyMulModGadget<FieldT> g(pb, in1[i][j][k], in1.max_value[i][j], eval_lc[j][k], eval_max_value[j],
                                           in1.modulus[j]);
                g.generate_r1cs_constraints();
                g.generate_r1cs_witness();
                out[i][j][k]            = g.out;
                out_max_values[i][j][k] = g.out_max_value;
            }
        }
    }
    for (size_t i = 0; i < out.size(); ++i) {
        out.max_value[i].resize(out[i].size());
        for (size_t j = 0; j < out[i].size(); ++j) {
            out.max_value[i][j] = get_max_field_element(out_max_values[i][j]);
        }
    }
    SetProofMetadata(ctxt_out, std::make_shared<LibsnarkProofMetadata>(out));
}

void LibsnarkProofSystem::ConstrainSquare(const Ciphertext<DCRTPoly>& ctxt, Ciphertext<DCRTPoly>& ctxt_out) {
    const auto in = *GetProofMetadata(ctxt);

    const auto num_limbs = in.modulus.size();
    assert(
        ctxt->GetElements()[0].GetNumOfElements() == ctxt_out->GetElements()[0].GetNumOfElements() &&
        "mismatch between number of limbs between ciphertext input and output. Are you using the FIXEDMANUAL scaling technique?");
    assert(in.size() == 2);

    LibsnarkProofMetadata tmp(1);
    for (size_t i = 0; i < tmp.size(); i++) {
        tmp[i] = vector<vector<pb_linear_combination<FieldT>>>(num_limbs);
    }
    tmp.modulus = in.modulus;
    constrain_mulmod_lazy(in, 0, in, 1, tmp, 0);

    LibsnarkProofMetadata out(3);
    for (size_t i = 0; i < out.size(); i++) {
        out[i] = vector<vector<pb_linear_combination<FieldT>>>(num_limbs);
    }
    out.modulus = in.modulus;
    constrain_mulmod_lazy(in, 0, in, 0, out, 0);
    constrain_addmod_lazy(tmp, 0, tmp, 0, out, 1);
    constrain_mulmod_lazy(in, 1, in, 1, out, 2);
    SetProofMetadata(ctxt_out, std::make_shared<LibsnarkProofMetadata>(out));
}

std::shared_ptr<LibsnarkProofMetadata> LibsnarkProofSystem::ConstrainPublicOutput(Ciphertext<DCRTPoly>& ciphertext) {
    const size_t num_poly  = ciphertext->GetElements().size();
    const size_t num_limbs = ciphertext->GetElements()[0].GetNumOfElements();
    LibsnarkProofMetadata out(num_poly);
    out.max_value = vector<vector<FieldT>>(num_poly);
    out.modulus   = vector<size_t>(num_limbs);

    for (size_t j = 0; j < num_limbs; j++) {
        out.modulus[j] = ciphertext->GetElements()[0].GetElementAtIndex(j).GetModulus().ConvertToInt<unsigned long>();
    }

    for (size_t i = 0; i < num_poly; i++) {
        const auto c_i   = ciphertext->GetElements()[i];
        out[i]           = vector<vector<pb_linear_combination<FieldT>>>(num_limbs);
        out.max_value[i] = vector<FieldT>(num_limbs);
        for (size_t j = 0; j < num_limbs; j++) {
            const auto c_ij     = c_i.GetElementAtIndex(j);
            const auto& v_ij    = c_ij.GetValues();
            out[i][j]           = vector<pb_linear_combination<FieldT>>(v_ij.GetLength());
            out.max_value[i][j] = FieldT(c_ij.GetModulus().ConvertToInt<size_t>()) - 1;
            for (size_t k = 0; k < v_ij.GetLength(); k++) {
                pb_variable<FieldT> tmp;
                tmp.allocate(pb, ciphertext->SerializedObjectName() + "[" + std::to_string(i) + "][" +
                                     std::to_string(j) + "][" + std::to_string(k) + "] (output)");
                out[i][j][k] = pb_linear_combination<FieldT>(tmp);
            }
        }
    }

    pb.set_input_sizes(pb.num_inputs() + (out.size() * out[0].size() * out[0][0].size()));
    return std::make_shared<LibsnarkProofMetadata>(out);
}

void LibsnarkProofSystem::FinalizeOutputConstraints(Ciphertext<DCRTPoly>& ctxt, const LibsnarkProofMetadata& out_vars) {
    // ctxt holds metadata for the output of the computation, vars holds the (public input) variables allocated at the beginning of the computation
    // We resolve all pending lazy mod-reductions, and add constraints binding vars to the output of the computation
    auto out           = *GetProofMetadata(ctxt);
    const auto modulus = out.modulus;

    assert(ctxt->GetElements().size() == out_vars.size());
    for (size_t i = 0; i < ctxt->GetElements().size(); ++i) {
        auto c_i = ctxt->GetElements()[i];
        for (size_t j = 0; j < c_i.GetNumOfElements(); ++j) {
            auto c_ij            = c_i.GetElementAtIndex(j);
            bool needs_reduction = gt_eq(out.max_value[i][j], FieldT(out.modulus[j]));
            vector<pb_variable<FieldT>> vars;
            vars.reserve(out[i][j].size());
            for (const auto& x : out_vars[i][j]) {
                assert(x.is_variable);
                vars.emplace_back(x.terms[0].index);
            }
            if (needs_reduction) {
                auto g = BatchGadget<FieldT, ModAssignGadget<FieldT>>(
                    pb, out[i][j], out.max_value[i][j], modulus[j], vars,
                    "finalize_output_constraints[" + std::to_string(i) + "][" + std::to_string(j) + "]");
                g.generate_r1cs_constraints();
                g.generate_r1cs_witness();
                out.max_value[i][j] = FieldT(out.modulus[j]) - 1;
            }
            else {
                for (size_t k = 0; k < c_ij.GetLength(); ++k) {
                    pb.add_r1cs_constraint(r1cs_constraint<FieldT>(out[i][j][k], 1, vars[k]),
                                           "finalize_output_constraints[" + std::to_string(i) + "][" +
                                               std::to_string(j) + "][" + std::to_string(k) + "]");
                }
            }
        }
    }
}

template <typename VecType>
void LibsnarkProofSystem::ConstrainSwitchModulus(
    const typename VecType::Integer& newModulus, const typename VecType::Integer& rootOfUnity,
    const typename VecType::Integer& modulusArb, const typename VecType::Integer& rootOfUnityArb,
    const PolyImpl<VecType>& in, const PolyImpl<VecType>& out, const vector<pb_linear_combination<FieldT>>& in_lc,
    const FieldT& in_max_value, vector<pb_linear_combination<FieldT>>& out_lc, FieldT& out_max_value) {
/**Switches the integers in the vector to values corresponding to the new
 * modulus.
 * Algorithm: Integer i, Old Modulus om, New Modulus nm,
 * delta = abs(om-nm):
 *  Case 1: om < nm
 *    if i > om/2
 *      i' = i + delta
 *  Case 2: om > nm
 *    if i > om/2
 *      i' = i - delta
 */
#ifdef PROOFSYSTEM_CHECK_STRICT
    for (int i = 0; i < out_lc.size(); i++) {
        in_lc[i].evaluate(pb);
        assert(lt_eq(pb.lc_val(in_lc[i]), in_max_value));
        assert(mod(pb.lc_val(in_lc[i]), FieldT(in.GetModulus().template ConvertToInt<unsigned long>())) ==
               FieldT(in[i].template ConvertToInt<unsigned long>()));
    }
#endif
    out_lc.resize(in_lc.size());

    auto oldModulus(in.GetModulus());
    auto oldModulusByTwo(oldModulus >> 1);
    auto diff((oldModulus > newModulus) ? (oldModulus - newModulus) : (newModulus - oldModulus));

    auto in_red_lc        = in_lc;
    auto in_red_max_value = in_max_value;
    auto old_mod_int      = oldModulus.template ConvertToInt<unsigned long>();
    if (gt_eq(in_max_value, FieldT(old_mod_int))) {
        for (usint i = 0; i < in.GetLength(); i++) {
            // We need to mod-reduce before continuing
            ModGadget<FieldT> g(pb, in_lc[i], in_max_value, old_mod_int);
            g.generate_r1cs_constraints();
            g.generate_r1cs_witness();
            in_red_lc[i] = g.out;
        }
        in_red_max_value = FieldT(old_mod_int) - 1;
    }
    if (newModulus > oldModulus) {
        for (usint i = 0; i < in.GetLength(); i++) {
            assert(oldModulusByTwo + diff < newModulus);
            // b == [ in <= oldModulusByTwo ]
            // out == b * in + (1-b) * (in + diff), which we simplify to out == in + (1-b) * diff
            is_less_than_constant_gadget<FieldT> g(pb, in_red_lc[i], in_red_max_value.as_bigint().num_bits(),
                                                   FieldT(oldModulusByTwo.template ConvertToInt<unsigned long>()) + 1);
            g.generate_r1cs_constraints();
            g.generate_r1cs_witness();
            out_lc[i].assign(pb,
                             in_red_lc[i] + (1 - g.less_or_eq) * FieldT(diff.template ConvertToInt<unsigned long>()));

            out_lc[i].evaluate(pb);
        }
    }
    else {
        assert(oldModulusByTwo > diff);  // If q/2 > q', then in-diff >= 0

        for (usint i = 0; i < in.GetLength(); i++) {
            // b == [ in <= oldModulusByTwo ]
            // tmp == b * in + (1-b) * (in - diff), which we simplify to tmp == in - (1-b) * diff
            // out == tmp (mod) newModulus
            is_less_than_constant_gadget<FieldT> g(pb, in_red_lc[i], in_red_max_value.as_bigint().num_bits(),
                                                   FieldT(oldModulusByTwo.template ConvertToInt<unsigned long>()) + 1);
            g.generate_r1cs_constraints();
            g.generate_r1cs_witness();
            pb_linear_combination<FieldT> tmp;
            tmp.assign(pb, in_red_lc[i] - (1 - g.less_or_eq) * FieldT(diff.template ConvertToInt<unsigned long>()));
            FieldT tmp_max_value =
                FieldT(std::max(oldModulusByTwo, oldModulus - 1 - diff).template ConvertToInt<unsigned long>());

            auto n        = in.GetValues()[i];
            auto sub_diff = (n > oldModulusByTwo) ? diff : 0;
            assert(n >= sub_diff);
            n = n.Sub(sub_diff);

            tmp.evaluate(pb);
#ifdef PROOFSYSTEM_CHECK_STRICT
//            assert(pb.lc_val(tmp) == FieldT(n.template ConvertToInt<unsigned long>()));
#endif
            ModGadget<FieldT> g_mod(pb, tmp, tmp_max_value, newModulus.template ConvertToInt<unsigned long>(), "",
                                    false);
            out_lc[i] = g_mod.out;
            g_mod.generate_r1cs_constraints();
            g_mod.generate_r1cs_witness();
        }
    }
    out_max_value = newModulus.template ConvertToInt<unsigned long>() - 1;
#ifdef PROOFSYSTEM_CHECK_STRICT
    for (int i = 0; i < out_lc.size(); i++) {
        out_lc[i].evaluate(pb);
        assert(lt_eq(pb.lc_val(out_lc[i]), out_max_value));
    }
#endif
}

template <typename VecType, typename VecType2>
void LibsnarkProofSystem::ConstrainNTT(const VecType& rootOfUnityTable, const VecType& preconRootOfUnityTable,
                                       const VecType2& element_in, const VecType2& element_out,
                                       const vector<pb_linear_combination<FieldT>>& in_lc, const FieldT& in_max_value,
                                       vector<pb_linear_combination<FieldT>>& out_lc, FieldT& out_max_value) {
    // Taken from NumberTheoreticTransformNat<VecType>::ForwardTransformToBitReverseInPlace
    auto element(element_in);
    assert(element.GetLength() == in_lc.size());
#ifdef PROOFSYSTEM_CHECK_STRICT
    for (int i = 0; i < in_lc.size(); i++) {
        in_lc[i].evaluate(pb);
        assert(lt_eq(pb.lc_val(in_lc[i]), in_max_value));
        assert(mod(pb.lc_val(in_lc[i]), FieldT(element.GetModulus().template ConvertToInt<unsigned long>())) ==
               FieldT(element[i].template ConvertToInt<unsigned long>()));
    }
#endif

    out_lc = vector<pb_linear_combination<FieldT>>(in_lc);
    vector<FieldT> out_max_values(out_lc.size(), in_max_value);

    assert(out_lc.size() == in_lc.size());
    using IntType = typename VecType::Integer;

    usint n         = element.GetLength();
    IntType modulus = element.GetModulus();
    unsigned long q = modulus.template ConvertToInt<unsigned long>();

    uint32_t indexOmega, indexHi;
    NativeInteger preconOmega;
    IntType omega, omegaFactor, loVal, hiVal, zero(0);

    usint t     = (n >> 1);
    usint logt1 = lbcrypto::GetMSB64(t);
    for (uint32_t m = 1; m < n; m <<= 1, t >>= 1, --logt1) {
        uint32_t j1, j2;
        for (uint32_t i = 0; i < m; ++i) {
            j1          = i << logt1;
            j2          = j1 + t;
            indexOmega  = m + i;
            omega       = rootOfUnityTable[indexOmega];
            preconOmega = preconRootOfUnityTable[indexOmega];
            for (uint32_t indexLo = j1; indexLo < j2; ++indexLo) {
                // omegaFactor = element[indexHi] * omega (mod q)
                // element_out[indexLo] = element[indexLo] + omegaFactor (mod q)
                // element_out[indexHi] = element[indexLo] - omegaFactor (mod q)
                indexHi     = indexLo + t;
                loVal       = element[indexLo];
                omegaFactor = element[indexHi];
                omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);

                auto oldLoVal = loVal;

                hiVal = loVal + omegaFactor;
                if (hiVal >= modulus) {
                    hiVal -= modulus;
                }

                if (loVal < omegaFactor) {
                    loVal += modulus;
                }
                loVal -= omegaFactor;

                // TODO: OPTIMIZEME: we might be able to use a LazyMulModGadget here in some cases, when out_max_values[indexHi] * omega is smaller than modulus
                assert(indexHi < out_lc.size());
                assert(indexHi < out_max_values.size());
                MulModGadget<FieldT> g1(pb, out_lc.at(indexHi), out_max_values.at(indexHi),
                                        FieldT(omega.template ConvertToInt<unsigned long>()), q);
                FieldT g1_out_max_value = FieldT(q) - 1;
                g1.generate_r1cs_constraints();
                g1.generate_r1cs_witness();

#ifdef PROOFSYSTEM_CHECK_STRICT
                PROOFSYSTEM_ASSERT_EQ(pb.val(g1.out), FieldT(omegaFactor.template ConvertToInt<unsigned long>()));
                out_lc[indexLo].evaluate(pb);
                out_lc[indexHi].evaluate(pb);
                PROOFSYSTEM_ASSERT_EQ(mod(pb.lc_val(out_lc[indexLo]), FieldT(q)),
                                      FieldT(element[indexLo].template ConvertToInt<unsigned long>()));
                PROOFSYSTEM_ASSERT_EQ(mod(pb.lc_val(out_lc[indexHi]), FieldT(q)),
                                      FieldT(element[indexHi].template ConvertToInt<unsigned long>()));

#endif
                LazyAddModGadget<FieldT> g2(pb, out_lc[indexLo], out_max_values[indexLo], g1.out, g1_out_max_value, q);
                //                LazyAddModGadget<FieldT> g2(pb, out_lc[indexLo], out_max_values[indexLo], FieldT(omegaFactor.template ConvertToInt<unsigned long>()), q);

                LazySubModGadget<FieldT> g3(pb, out_lc[indexLo], out_max_values[indexLo], g1.out, g1_out_max_value, q);

                g2.generate_r1cs_constraints();
                g2.generate_r1cs_witness();
                g3.generate_r1cs_constraints();
                g3.generate_r1cs_witness();

                out_lc[indexLo]         = g2.out;
                out_max_values[indexLo] = g2.out_max_value;
                out_lc[indexHi]         = g3.out;
                out_max_values[indexHi] = g3.out_max_value;

                element[indexLo] = hiVal;
                element[indexHi] = loVal;

                assert(element[indexLo] == (oldLoVal + omegaFactor).Mod(q));
                assert(element[indexHi] == (oldLoVal + q - omegaFactor).Mod(q));

#ifdef PROOFSYSTEM_CHECK_STRICT
                out_lc[indexLo].evaluate(pb);
                out_lc[indexHi].evaluate(pb);
                assert(lt_eq(pb.lc_val(out_lc[indexLo]), out_max_values[indexLo]));
                assert(lt_eq(pb.lc_val(out_lc[indexHi]), out_max_values[indexHi]));

                PROOFSYSTEM_ASSERT_EQ(mod(pb.lc_val(out_lc[indexLo]), FieldT(q)),
                                      FieldT(element[indexLo].template ConvertToInt<unsigned long>()));
                PROOFSYSTEM_ASSERT_EQ(mod(pb.lc_val(out_lc[indexHi]), FieldT(q)),
                                      FieldT(element[indexHi].template ConvertToInt<unsigned long>()));
#endif
            }
        }
    }

    // Set out_max_value to max of all out_max_values for soundness
    out_max_value = 0;
    for (size_t i = 0; i < n; i++) {
        if (lt(out_max_value, out_max_values[i])) {
            out_max_value = out_max_values[i];
        }
    }

#ifdef PROOFSYSTEM_CHECK_STRICT
    for (int i = 0; i < n; i++) {
        assert(element[i] == element_out[i]);
    }
    for (int i = 0; i < out_lc.size(); i++) {
        out_lc[i].evaluate(pb);
        assert(lt_eq(pb.lc_val(out_lc[i]), out_max_value));
        assert(mod(pb.lc_val(out_lc[i]), FieldT(q)) == FieldT(element_out[i].template ConvertToInt<unsigned long>()));
    }
#endif
}

template <typename IntType, typename VecType, typename VecType2>
void LibsnarkProofSystem::ConstrainINTT(const VecType& rootOfUnityInverseTable,
                                        const VecType& preconRootOfUnityInverseTable, const IntType& cycloOrderInv,
                                        const IntType& preconCycloOrderInv, const VecType2& element_in,
                                        const VecType2& element_out, const vector<pb_linear_combination<FieldT>>& in_lc,
                                        const FieldT& in_max_value, vector<pb_linear_combination<FieldT>>& out_lc,
                                        FieldT& out_max_value) {
    // Taken from NumberTheoreticTransformNat<VecType>::InverseTransformFromBitReverseInPlace

    auto element(element_in);
#ifdef PROOFSYSTEM_CHECK_STRICT
    for (int i = 0; i < in_lc.size(); i++) {
        in_lc[i].evaluate(pb);
        assert(lt_eq(pb.lc_val(in_lc[i]), in_max_value));
        assert(mod(pb.lc_val(in_lc[i]), FieldT(element.GetModulus().template ConvertToInt<unsigned long>())) ==
               FieldT(element[i].template ConvertToInt<unsigned long>()));
    }
#endif

    out_lc = vector<pb_linear_combination<FieldT>>(in_lc);
    vector<FieldT> out_max_values(out_lc.size(), in_max_value);

    usint n = element.GetLength();

    IntType modulus = element.GetModulus();
    unsigned long q = modulus.template ConvertToInt<unsigned long>();

    IntType loVal, hiVal, omega, omegaFactor;
    NativeInteger preconOmega;
    usint i, m, j1, j2, indexOmega, indexLo, indexHi;

    usint t     = 1;
    usint logt1 = 1;
    for (m = (n >> 1); m >= 1; m >>= 1) {
        for (i = 0; i < m; ++i) {
            j1          = i << logt1;
            j2          = j1 + t;
            indexOmega  = m + i;
            omega       = rootOfUnityInverseTable[indexOmega];
            preconOmega = preconRootOfUnityInverseTable[indexOmega];

            for (indexLo = j1; indexLo < j2; ++indexLo) {
                // element_out[indexLo] = element[indexLo] + element[indexHi] (mod q)
                // element_out[indexHi] = (element[indexLo] - element[indexHi]) * omega (mod q)
                indexHi = indexLo + t;

                hiVal = element[indexHi];
                loVal = element[indexLo];

                omegaFactor = loVal;
                if (omegaFactor < hiVal) {
                    omegaFactor += modulus;
                }

                omegaFactor -= hiVal;

                loVal += hiVal;
                if (loVal >= modulus) {
                    loVal -= modulus;
                }

                omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);

                LazyAddModGadget<FieldT> g1(this->pb, out_lc[indexLo], out_max_values[indexLo], out_lc[indexHi],
                                            out_max_values[indexHi], q);

                // If out_lc[indexHi] is > modulus, we cannot use the lazy sub-mod gadget as-is, and we need to reduce out_lc[indexHi] mod modulus first
                pb_linear_combination<FieldT> hi_reduced;
                FieldT hi_reduced_max_value;
                if (lt_eq(out_max_values[indexHi], FieldT(q))) {
                    hi_reduced           = out_lc[indexHi];
                    hi_reduced_max_value = out_max_values[indexHi];
                }
                else {
                    ModGadget<FieldT> g_mod(this->pb, out_lc[indexHi], out_max_values[indexHi], q);
                    g_mod.generate_r1cs_constraints();
                    g_mod.generate_r1cs_witness();
                    hi_reduced           = g_mod.out;
                    hi_reduced_max_value = FieldT(q) - 1;
                }

                LazySubModGadget<FieldT> g2(this->pb, out_lc[indexLo], out_max_values[indexLo], hi_reduced,
                                            hi_reduced_max_value, q);
                LazyMulModGadget<FieldT> g3(this->pb, g2.out, g2.out_max_value, FieldT(omega.ConvertToInt()), q);

                g1.generate_r1cs_constraints();
                g1.generate_r1cs_witness();

                g2.generate_r1cs_constraints();
                g2.generate_r1cs_witness();
                g3.generate_r1cs_constraints();
                g3.generate_r1cs_witness();

                out_lc[indexLo]         = g1.out;
                out_max_values[indexLo] = g1.out_max_value;
                //                out_lc[indexLo].evaluate(pb);

                out_lc[indexHi]         = g3.out;
                out_max_values[indexHi] = g3.out_max_value;
                //                out_lc[indexHi].evaluate(pb);

                element[indexLo] = loVal;
                element[indexHi] = omegaFactor;

#ifdef PROOFSYSTEM_CHECK_STRICT
                out_lc[indexLo].evaluate(pb);
                out_lc[indexHi].evaluate(pb);
                assert(lt_eq(pb.lc_val(out_lc[indexLo]), out_max_values[indexLo]));
                assert(lt_eq(pb.lc_val(out_lc[indexHi]), out_max_values[indexHi]));

                PROOFSYSTEM_ASSERT_EQ(mod(pb.lc_val(out_lc[indexLo]), FieldT(q)),
                                      FieldT(element[indexLo].template ConvertToInt<unsigned long>()));
                PROOFSYSTEM_ASSERT_EQ(mod(pb.lc_val(out_lc[indexHi]), FieldT(q)),
                                      FieldT(element[indexHi].template ConvertToInt<unsigned long>()));
#endif
            }
        }
        t <<= 1;
        logt1++;
    }

    out_max_value = 0;
    for (i = 0; i < n; i++) {
        element[i].ModMulFastConstEq(cycloOrderInv, modulus, preconCycloOrderInv);
        LazyMulModGadget<FieldT> g(this->pb, out_lc[i], out_max_values[i], FieldT(cycloOrderInv.ConvertToInt()), q);
        g.generate_r1cs_constraints();
        g.generate_r1cs_witness();
        out_lc[i]         = g.out;
        out_max_values[i] = g.out_max_value;

        // Set out_max_value to max of all out_max_values for soundness
        if (lt(out_max_value, out_max_values[i])) {
            out_max_value = out_max_values[i];
        }
    }

#ifdef PROOFSYSTEM_CHECK_STRICT
    for (int i = 0; i < n; i++) {
        assert(element[i] == element_out[i]);
    }
    for (int i = 0; i < out_lc.size(); i++) {
        out_lc[i].evaluate(pb);
        assert(lt_eq(pb.lc_val(out_lc[i]), out_max_value));
        assert(mod(pb.lc_val(out_lc[i]), FieldT(q)) == FieldT(element_out[i].template ConvertToInt<unsigned long>()));
    }
#endif
}

template <typename VecType2>
void LibsnarkProofSystem::ConstrainSetFormat(const Format format, const VecType2& in, const VecType2& out,
                                             const vector<pb_linear_combination<FieldT>>& in_lc,
                                             const FieldT& in_max_value, vector<pb_linear_combination<FieldT>>& out_lc,
                                             FieldT& out_max_value) {
    assert(out.GetFormat() == format);
    assert(in.GetLength() == out.GetLength());
    assert(in.GetLength() == in_lc.size());
    if (in.GetFormat() == format) {
        assert(in == out);
        out_lc        = in_lc;
        out_max_value = in_max_value;
        return;
    }
    assert(in.GetFormat() != format);

    using namespace intnat;

    auto rootOfUnity    = out.GetRootOfUnity();
    auto CycloOrder     = out.GetCyclotomicOrder();
    const auto& modulus = in.GetModulus();
    usint CycloOrderHf  = (CycloOrder >> 1);

    using VecType = typename DCRTPoly::PolyType::Vector;

    if (format == COEFFICIENT) {
        auto mapSearch = ChineseRemainderTransformFTTNat<VecType>::m_rootOfUnityReverseTableByModulus.find(modulus);
        if (mapSearch == ChineseRemainderTransformFTTNat<VecType>::m_rootOfUnityReverseTableByModulus.end() ||
            mapSearch->second.GetLength() != CycloOrderHf) {
            ChineseRemainderTransformFTTNat<VecType>().PreCompute(rootOfUnity, CycloOrder, modulus);
        }
        usint msb = lbcrypto::GetMSB64(CycloOrderHf - 1);

        ConstrainINTT(
            ChineseRemainderTransformFTTNat<VecType>::m_rootOfUnityInverseReverseTableByModulus[modulus],
            ChineseRemainderTransformFTTNat<VecType>::m_rootOfUnityInversePreconReverseTableByModulus[modulus],
            ChineseRemainderTransformFTTNat<VecType>::m_cycloOrderInverseTableByModulus[modulus][msb],
            ChineseRemainderTransformFTTNat<VecType>::m_cycloOrderInversePreconTableByModulus[modulus][msb], in, out,
            in_lc, in_max_value, out_lc, out_max_value);
    }
    else {
        auto mapSearch = ChineseRemainderTransformFTTNat<VecType>::m_rootOfUnityReverseTableByModulus.find(modulus);
        if (mapSearch == ChineseRemainderTransformFTTNat<VecType>::m_rootOfUnityReverseTableByModulus.end() ||
            mapSearch->second.GetLength() != CycloOrderHf) {
            ChineseRemainderTransformFTTNat<VecType>().PreCompute(rootOfUnity, CycloOrder, modulus);
        }

        ConstrainNTT(ChineseRemainderTransformFTTNat<VecType>::m_rootOfUnityReverseTableByModulus[modulus],
                     ChineseRemainderTransformFTTNat<VecType>::m_rootOfUnityPreconReverseTableByModulus[modulus], in,
                     out, in_lc, in_max_value, out_lc, out_max_value);
    }
}

void LibsnarkProofSystem::ConstrainSetFormat(const Format format, const DCRTPoly& in, const DCRTPoly& out,
                                             const vector<vector<pb_linear_combination<FieldT>>>& in_lc,
                                             const vector<FieldT>& in_max_value,
                                             vector<vector<pb_linear_combination<FieldT>>>& out_lc,
                                             vector<FieldT>& out_max_value) {
    assert(out.GetFormat() == format);
    size_t n = in.GetNumOfElements();
    assert(out.GetNumOfElements() == n);
    assert(in_lc.size() == n);
    assert(in_max_value.size() == n);
    assert(in.GetFormat() != format);
    if (in.GetFormat() == format) {
        assert(in == out);
        out_lc        = in_lc;
        out_max_value = in_max_value;
        return;
    }
    assert(in.GetFormat() != format);
    out_lc.resize(n);
    out_max_value.resize(n);

    for (size_t i = 0; i < n; i++) {
        ConstrainSetFormat(format, in.GetElementAtIndex(i), out.GetElementAtIndex(i), in_lc[i], in_max_value[i],
                           out_lc[i], out_max_value[i]);
    }
}
template <typename DCRTPoly>
void LibsnarkProofSystem::ConstrainRescale(const Ciphertext<DCRTPoly>& ctxt_in, Ciphertext<DCRTPoly>& ctxt_out) {
    auto in                 = *GetProofMetadata(ctxt_in);
    const size_t num_polys  = in.size();
    const size_t num_levels = in[0].size();
    assert(ctxt_in->GetElements().size() == num_polys);
    assert(ctxt_in->GetElements()[0].GetNumOfElements() == num_levels);

    LibsnarkProofMetadata out(num_polys);
    out.modulus.resize(num_levels - 1);
    for (size_t j = 0; j < num_levels - 1; ++j) {
        out.modulus[j] = in.modulus[j];
    }
    for (size_t i = 0; i < num_polys; ++i) {
        out[i].resize(num_levels - 1);
        out.max_value[i].resize(num_levels - 1);
        for (size_t j = 0; j < num_levels - 1; ++j) {
            out[i][j]           = in[i][j];
            out.max_value[i][j] = in.max_value[i][j];
        }

        vector<pb_linear_combination<FieldT>> last_lc;
        FieldT last_max_value;
        auto lastPoly = ctxt_in->GetElements()[i].GetElementAtIndex(num_levels - 1);
        auto lastPolyCoef(lastPoly);
        lastPolyCoef.SetFormat(Format::COEFFICIENT);

        ConstrainSetFormat(Format::COEFFICIENT, lastPoly, lastPolyCoef, in[i][num_levels - 1],
                           in.max_value[i][num_levels - 1], last_lc, last_max_value);

        DCRTPoly extra(ctxt_in->GetElements()[i].GetParams(), COEFFICIENT, true);
        vector<vector<pb_linear_combination<FieldT>>> extra_lc(extra.GetNumOfElements());
        vector<FieldT> extra_max_value(extra.GetNumOfElements());
        for (usint j = 0; j < extra.GetNumOfElements(); j++) {
            auto temp         = lastPoly;
            const auto newMod = ctxt_in->GetElements()[0].GetElementAtIndex(j).GetModulus();
            const auto newRoU = ctxt_in->GetElements()[0].GetElementAtIndex(j).GetRootOfUnity();

            temp.SwitchModulus(newMod, newRoU, 0, 0);
            ConstrainSwitchModulus(newMod, newRoU, 0, 0, lastPoly, temp, last_lc, last_max_value, extra_lc[j],
                                   extra_max_value[j]);
            //            extra.m_vectors[j] = (temp *= QlQlInvModqlDivqlModq[j]);
        }

        //        if (this->GetFormat() == Format::EVALUATION)
        //            extra.SetFormat(Format::EVALUATION);
        auto extra_eval(extra);
        extra_eval.SetFormat(Format::EVALUATION);
        vector<vector<pb_linear_combination<FieldT>>> extra_eval_lc(extra.GetNumOfElements());
        vector<FieldT> extra_eval_max_value(extra.GetNumOfElements());
        ConstrainSetFormat(Format::EVALUATION, extra, extra_eval, extra_lc, extra_max_value, extra_eval_lc,
                           extra_eval_max_value);

        //        for (usint i = 0; i < m_vectors.size(); i++) {
        //            m_vectors[i] *= qlInvModq[i];
        //            m_vectors[i] += extra.m_vectors[i];
        //        }

        for (size_t j = 0; j < out[i].size(); j++) {
            FieldT curr_max = 0;
            for (size_t k = 0; k < out[i][j].size(); k++) {
                LazyAddModGadget<FieldT> g(pb, out[i][j][k], out.max_value.at(i).at(j), extra_eval_lc.at(j).at(k),
                                           extra_eval_max_value.at(j), out.modulus.at(j));
                g.generate_r1cs_constraints();
                g.generate_r1cs_witness();
                out[i][j][k] = g.out;
                if (lt(curr_max, g.out_max_value)) {
                    curr_max = g.out_max_value;
                }
            }
            out.max_value[i][j] = curr_max;
        }

        //        this->SetFormat(Format::EVALUATION);
        assert(ctxt_out->GetElements()[i].GetFormat() == Format::EVALUATION);
        //        auto curr_eval(ctxt_out->GetElements()[i]);
        //        curr_eval.SetFormat(Format::EVALUATION);
        //        ConstrainSetFormat(Format::EVALUATION, ctxt_out->GetElements()[i], curr_eval, out[i], out.max_value[i], out[i],
        //                           out.max_value[i]);
    }

    SetProofMetadata(ctxt_out, std::make_shared<LibsnarkProofMetadata>(out));
}

template <typename DCRTPoly>
void LibsnarkProofSystem::ConstrainKeySwitchPrecomputeCore(
    const DCRTPoly& in, const std::shared_ptr<CryptoParametersBase<DCRTPoly>>& cryptoParamsBase,
    const std::shared_ptr<std::vector<DCRTPoly>>& out, const vector<vector<pb_linear_combination<FieldT>>>& in_lc,
    const vector<FieldT>& in_max_value, vector<vector<vector<pb_linear_combination<FieldT>>>>& out_lc,
    vector<vector<FieldT>>& out_max_value) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersRNS>(cryptoParamsBase);

    const size_t num_levels = in_lc.size();

    out_lc.resize(num_levels);
    out_max_value.resize(num_levels);
    for (size_t i = 0; i < num_levels; ++i) {
        out_lc[i].resize(num_levels);
        out_max_value[i].resize(num_levels);
    }

    // Taken from DCRTPolyImpl<VecType>::CRTDecompose

    //    auto decomposed = c.CRTDecompose(digitSize);
    //    return std::make_shared<std::vector<DCRTPoly>>(decomposed.begin(), decomposed.end());

    uint32_t digitSize = cryptoParams->GetDigitSize();

    // used to store the number of digits for each small modulus
    uint32_t nWindows = 0;

    std::vector<usint> arrWindows;
    auto baseBits = digitSize;
    assert(baseBits == 0);

    /*     if (baseBits > 0) {
            nWindows = 0;

            // creates an array of digits up to a certain tower
            for (usint i = 0; i < m_vectors.size(); i++) {
                usint nBits      = m_vectors[i].GetModulus().GetLengthForBase(2);
                usint curWindows = nBits / baseBits;
                if (nBits % baseBits > 0)
                    curWindows++;
                arrWindows.push_back(nWindows);
                nWindows += curWindows;
            }
        }
        else {*/
    nWindows = in.GetNumOfElements();
    //    }

    using DCRTPolyType = DCRTPoly;
    std::vector<DCRTPolyType> result(nWindows);
    result = *out;

    DCRTPolyType input = in.Clone();
    input.SetFormat(Format::COEFFICIENT);
    vector<vector<pb_linear_combination<FieldT>>> in_coeff_lc;
    vector<FieldT> in_coeff_max_value;
    ConstrainSetFormat(Format::COEFFICIENT, in, input, in_lc, in_max_value, in_coeff_lc, in_coeff_max_value);

    // out[k] holds a representation of the k-th limb of in, i.e., out[k] = f(in[k])
    for (usint i = 0; i < num_levels; i++) {
        if (baseBits == 0) {
            //            DCRTPolyType currentDCRTPoly = input.Clone();

            for (usint k = 0; k < num_levels; k++) {
                auto temp(input.GetElementAtIndex(i));
                auto old_temp(temp);
                auto old_temp_lc        = in_coeff_lc[i];
                auto old_temp_max_value = in_coeff_max_value[i];
                vector<pb_linear_combination<FieldT>> temp_lc;
                FieldT temp_max_value;
                if (i != k) {
                    temp.SwitchModulus(input.GetElementAtIndex(k).GetModulus(),
                                       input.GetElementAtIndex(k).GetRootOfUnity(), 0, 0);
                    ConstrainSwitchModulus(input.GetElementAtIndex(k).GetModulus(),
                                           input.GetElementAtIndex(k).GetRootOfUnity(), 0, 0, old_temp, temp,
                                           old_temp_lc, old_temp_max_value, temp_lc, temp_max_value);

                    // temp.SetFormat(Format::EVALUATION);
                    ConstrainSetFormat(Format::EVALUATION, temp, (*out)[i].GetElementAtIndex(k), temp_lc,
                                       temp_max_value, out_lc[i][k], out_max_value[i][k]);
                }
                else {  // saves an extra NTT
                    // currentDCRTPoly.m_vectors[k] = this->m_vectors[k];
                    // auto curr_coef = input.GetElementAtIndex(k);
                    // auto curr_eval(curr_coef);
                    // curr_eval.SetFormat(Format::EVALUATION);

                    // currentDCRTPoly.m_vectors[k].SetFormat(Format::EVALUATION);
                    ConstrainSetFormat(Format::EVALUATION, input.GetElementAtIndex(k), (*out)[i].GetElementAtIndex(k),
                                       in_coeff_lc[k], in_coeff_max_value[k], out_lc[i][k], out_max_value[i][k]);
                    assert((*out)[i].GetElementAtIndex(k) == in.GetElementAtIndex(k));
                }
            }
            // currentDCRTPoly.m_format = Format::EVALUATION;
            // result[i] = std::move(currentDCRTPoly);
        }
    }
}

template <typename DCRTPoly>
void LibsnarkProofSystem::ConstrainFastKeySwitchCore(const EvalKey<DCRTPoly>& evalKey,
                                                     const std::shared_ptr<typename DCRTPoly::Params>& paramsQl,
                                                     const vector<vector<vector<pb_linear_combination<FieldT>>>>& in_lc,
                                                     const vector<vector<FieldT>>& in_max_value,
                                                     vector<vector<vector<pb_linear_combination<FieldT>>>>& out_lc,
                                                     vector<vector<FieldT>>& out_max_value) {
    const auto n = in_lc[0].size();

    std::vector<DCRTPoly> bv(evalKey->GetBVector());
    std::vector<DCRTPoly> av(evalKey->GetAVector());

    auto sizeQ    = bv[0].GetParams()->GetParams().size();
    auto sizeQl   = paramsQl->GetParams().size();
    size_t diffQl = sizeQ - sizeQl;

    for (size_t k = 0; k < bv.size(); k++) {
        av[k].DropLastElements(diffQl);
        bv[k].DropLastElements(diffQl);
    }

    out_lc        = vector<vector<vector<pb_linear_combination<FieldT>>>>(in_lc);
    out_max_value = (in_max_value);

    // av, bv are public constants, digits are private variables
    //    DCRTPoly ct1 = (av[0] * (*digits)[0]);
    //    DCRTPoly ct0 = (bv[0] * (*digits)[0]);

    out_lc = vector<vector<vector<pb_linear_combination<FieldT>>>>(2);
    out_lc[0].resize(n);
    out_lc[1].resize(n);

    out_max_value = vector<vector<FieldT>>(2);
    out_max_value[0].resize(n);
    out_max_value[1].resize(n);

    auto ct_max_value = vector<vector<vector<FieldT>>>(2);
    ct_max_value[0].resize(n);
    ct_max_value[1].resize(n);

    for (size_t j = 0; j < n; ++j) {
        out_lc[0][j].resize(in_lc[0][j].size());
        out_lc[1][j].resize(in_lc[0][j].size());
        ct_max_value[0][j].resize(in_lc[0][j].size());
        ct_max_value[1][j].resize(in_lc[0][j].size());
        size_t modulus = av[0].GetElementAtIndex(j).GetModulus().ConvertToInt();
        for (size_t k = 0; k < in_lc[0][j].size(); ++k) {
            auto av_0jk = av[0].GetElementAtIndex(j).GetValues()[k];
            LazyMulModGadget<FieldT> g1(pb, in_lc[0][j][k], in_max_value[0][j], FieldT(av_0jk.ConvertToInt()), modulus);
            g1.generate_r1cs_constraints();
            g1.generate_r1cs_witness();
            out_lc[1][j][k]       = g1.out;
            ct_max_value[1][j][k] = g1.out_max_value;

            auto bv_0jk = bv[0].GetElementAtIndex(j).GetValues()[k];
            LazyMulModGadget<FieldT> g0(pb, in_lc[0][j][k], in_max_value[0][j], FieldT(bv_0jk.ConvertToInt()), modulus);
            g0.generate_r1cs_constraints();
            g0.generate_r1cs_witness();
            out_lc[0][j][k]       = g0.out;
            ct_max_value[0][j][k] = g0.out_max_value;
        }
    }

    for (usint i = 1; i < in_lc.size(); ++i) {
        //        ct0 += (bv[i] * (*digits)[i]);
        //        ct1 += (av[i] * (*digits)[i]);

        for (size_t j = 0; j < n; ++j) {
            out_lc[0][j].resize(in_lc[0][j].size());
            out_lc[1][j].resize(in_lc[0][j].size());
            size_t modulus = av[0].GetElementAtIndex(j).GetModulus().ConvertToInt();
            for (size_t k = 0; k < in_lc[0][j].size(); ++k) {
                auto bv_ijk = bv[i].GetElementAtIndex(j).GetValues()[k];
                LazyMulModGadget<FieldT> g0(pb, in_lc[i][j][k], in_max_value[i][j], FieldT(bv_ijk.ConvertToInt()),
                                            modulus);
                g0.generate_r1cs_constraints();
                g0.generate_r1cs_witness();
                LazyAddModGadget<FieldT> g0_add(pb, g0.out, g0.out_max_value, out_lc[0][j][k], ct_max_value[0][j][k],
                                                modulus);
                g0_add.generate_r1cs_constraints();
                g0_add.generate_r1cs_witness();
                out_lc[0][j][k]       = g0_add.out;
                ct_max_value[0][j][k] = g0_add.out_max_value;

                auto av_ijk = av[i].GetElementAtIndex(j).GetValues()[k];
                LazyMulModGadget<FieldT> g1(pb, in_lc[i][j][k], in_max_value[i][j], FieldT(av_ijk.ConvertToInt()),
                                            modulus);
                g1.generate_r1cs_constraints();
                g1.generate_r1cs_witness();
                LazyAddModGadget<FieldT> g1_add(pb, g1.out, g1.out_max_value, out_lc[1][j][k], ct_max_value[1][j][k],
                                                modulus);
                g0_add.generate_r1cs_constraints();
                g0_add.generate_r1cs_witness();
                out_lc[1][j][k]       = g1_add.out;
                ct_max_value[1][j][k] = g1_add.out_max_value;
            }
        }
    }
    for (size_t i = 0; i < ct_max_value.size(); ++i) {
        for (size_t j = 0; j < ct_max_value[i].size(); ++j) {
            out_max_value[i][j] = 0;
            for (size_t k = 0; k < ct_max_value[i][j].size(); ++k) {
                if (gt(ct_max_value[i][j][k], out_max_value[i][j])) {
                    out_max_value[i][j] = ct_max_value[i][j][k];
                }
            }
        }
    }
}

template <typename DCRTPoly>
void LibsnarkProofSystem::ConstrainFastKeySwitchCore(
    const std::shared_ptr<std::vector<DCRTPoly>>& digits, const EvalKey<DCRTPoly>& evalKey,
    const std::shared_ptr<typename DCRTPoly::Params>& paramsQl, std::shared_ptr<std::vector<DCRTPoly>>& out,
    const vector<vector<vector<pb_linear_combination<FieldT>>>>& in_lc, const vector<vector<FieldT>>& in_max_value,
    vector<vector<vector<pb_linear_combination<FieldT>>>>& out_lc, vector<vector<FieldT>>& out_max_value) {
    std::vector<DCRTPoly> bv(evalKey->GetBVector());
    std::vector<DCRTPoly> av(evalKey->GetAVector());

    auto sizeQ    = bv[0].GetParams()->GetParams().size();
    auto sizeQl   = paramsQl->GetParams().size();
    size_t diffQl = sizeQ - sizeQl;

    for (size_t k = 0; k < bv.size(); k++) {
        av[k].DropLastElements(diffQl);
        bv[k].DropLastElements(diffQl);
    }

    out_lc        = vector<vector<vector<pb_linear_combination<FieldT>>>>(in_lc);
    out_max_value = (in_max_value);

    // av, bv are public constants, digits are private variables
    DCRTPoly ct1 = (av[0] * (*digits)[0]);
    DCRTPoly ct0 = (bv[0] * (*digits)[0]);
    //    DCRTPoly ct1((*digits)[0]);
    //    DCRTPoly ct0((*digits)[0]);
    //    ct1 *= av[0];
    //    ct0 *= bv[0];

    out_lc = vector<vector<vector<pb_linear_combination<FieldT>>>>(2);
    out_lc[0].resize((*digits)[0].GetNumOfElements());
    out_lc[1].resize((*digits)[0].GetNumOfElements());

    out_max_value = vector<vector<FieldT>>(2);
    out_max_value[0].resize((*digits)[0].GetNumOfElements());
    out_max_value[1].resize((*digits)[0].GetNumOfElements());

    auto ct_max_value = vector<vector<vector<FieldT>>>(2);
    ct_max_value[0].resize((*digits)[0].GetNumOfElements());
    ct_max_value[1].resize((*digits)[0].GetNumOfElements());

    for (size_t j = 0; j < (*digits)[0].GetNumOfElements(); ++j) {
        out_lc[0][j].resize((*digits)[0].GetElementAtIndex(j).GetLength());
        out_lc[1][j].resize((*digits)[0].GetElementAtIndex(j).GetLength());
        ct_max_value[0][j].resize((*digits)[0].GetElementAtIndex(j).GetLength());
        ct_max_value[1][j].resize((*digits)[0].GetElementAtIndex(j).GetLength());
        size_t modulus = (*digits)[0].GetElementAtIndex(j).GetModulus().ConvertToInt();
        for (size_t k = 0; k < (*digits)[0].GetElementAtIndex(j).GetLength(); ++k) {
#ifdef PROOFSYSTEM_CHECK_STRICT
            in_lc[0][j][k].evaluate(pb);
            assert(lt_eq(pb.lc_val(in_lc[0][j][k]), in_max_value[0][j]));
            PROOFSYSTEM_ASSERT_EQ(mod(pb.lc_val(in_lc[0][j][k]), FieldT(modulus)),
                                  (*digits)[0].GetElementAtIndex(j).GetValues()[k].ConvertToInt());
#endif
            auto av_0jk = av[0].GetElementAtIndex(j).GetValues()[k];
            LazyMulModGadget<FieldT> g1(pb, in_lc[0][j][k], in_max_value[0][j], FieldT(av_0jk.ConvertToInt()), modulus);
            g1.generate_r1cs_constraints();
            g1.generate_r1cs_witness();
            out_lc[1][j][k]       = g1.out;
            ct_max_value[1][j][k] = g1.out_max_value;
#ifdef PROOFSYSTEM_CHECK_STRICT
            out_lc[1][j][k].evaluate(pb);
            assert(lt_eq(pb.lc_val(out_lc[1][j][k]), ct_max_value[1][j][k]));
            PROOFSYSTEM_ASSERT_EQ(mod(pb.lc_val(out_lc[1][j][k]), FieldT(modulus)),
                                  ct1.GetElementAtIndex(j).GetValues()[k].ConvertToInt());
#endif

            auto bv_0jk = bv[0].GetElementAtIndex(j).GetValues()[k];
            LazyMulModGadget<FieldT> g0(pb, in_lc[0][j][k], in_max_value[0][j], FieldT(bv_0jk.ConvertToInt()), modulus);
            g0.generate_r1cs_constraints();
            g0.generate_r1cs_witness();
            out_lc[0][j][k]       = g0.out;
            ct_max_value[0][j][k] = g0.out_max_value;
#ifdef PROOFSYSTEM_CHECK_STRICT
            out_lc[0][j][k].evaluate(pb);
            assert(lt_eq(pb.lc_val(out_lc[0][j][k]), ct_max_value[0][j][k]));
            PROOFSYSTEM_ASSERT_EQ(mod(pb.lc_val(out_lc[0][j][k]), FieldT(modulus)),
                                  ct0.GetElementAtIndex(j).GetValues()[k].ConvertToInt());
#endif
        }
    }

    for (usint i = 1; i < (*digits).size(); ++i) {
        ct0 += (bv[i] * (*digits)[i]);
        ct1 += (av[i] * (*digits)[i]);

        for (size_t j = 0; j < (*digits)[0].GetNumOfElements(); ++j) {
            out_lc[0][j].resize((*digits)[0].GetElementAtIndex(j).GetLength());
            out_lc[1][j].resize((*digits)[0].GetElementAtIndex(j).GetLength());
            size_t modulus = (*digits)[0].GetElementAtIndex(j).GetModulus().ConvertToInt();
            for (size_t k = 0; k < (*digits)[0].GetElementAtIndex(j).GetLength(); ++k) {
                auto bv_ijk = bv[i].GetElementAtIndex(j).GetValues()[k];
                LazyMulModGadget<FieldT> g0(pb, in_lc[i][j][k], in_max_value[i][j], FieldT(bv_ijk.ConvertToInt()),
                                            modulus);
                g0.generate_r1cs_constraints();
                g0.generate_r1cs_witness();
                LazyAddModGadget<FieldT> g0_add(pb, g0.out, g0.out_max_value, out_lc[0][j][k], ct_max_value[0][j][k],
                                                modulus);
                g0_add.generate_r1cs_constraints();
                g0_add.generate_r1cs_witness();
                out_lc[0][j][k]       = g0_add.out;
                ct_max_value[0][j][k] = g0_add.out_max_value;
#ifdef PROOFSYSTEM_CHECK_STRICT
                out_lc[0][j][k].evaluate(pb);
                assert(lt_eq(pb.lc_val(out_lc[0][j][k]), ct_max_value[0][j][k]));
                PROOFSYSTEM_ASSERT_EQ(mod(pb.lc_val(out_lc[0][j][k]), FieldT(modulus)),
                                      ct0.GetElementAtIndex(j).GetValues()[k].ConvertToInt());
#endif

                auto av_ijk = av[i].GetElementAtIndex(j).GetValues()[k];
                LazyMulModGadget<FieldT> g1(pb, in_lc[i][j][k], in_max_value[i][j], FieldT(av_ijk.ConvertToInt()),
                                            modulus);
                g1.generate_r1cs_constraints();
                g1.generate_r1cs_witness();
                LazyAddModGadget<FieldT> g1_add(pb, g1.out, g1.out_max_value, out_lc[1][j][k], ct_max_value[1][j][k],
                                                modulus);
                g0_add.generate_r1cs_constraints();
                g0_add.generate_r1cs_witness();
                out_lc[1][j][k]       = g1_add.out;
                ct_max_value[1][j][k] = g1_add.out_max_value;
#ifdef PROOFSYSTEM_CHECK_STRICT
                out_lc[1][j][k].evaluate(pb);
                assert(lt_eq(pb.lc_val(out_lc[1][j][k]), ct_max_value[1][j][k]));
                PROOFSYSTEM_ASSERT_EQ(mod(pb.lc_val(out_lc[1][j][k]), FieldT(modulus)),
                                      ct1.GetElementAtIndex(j).GetValues()[k].ConvertToInt());
#endif
            }
        }
    }
    for (size_t i = 0; i < ct_max_value.size(); ++i) {
        for (size_t j = 0; j < ct_max_value[i].size(); ++j) {
            out_max_value[i][j] = 0;
            for (size_t k = 0; k < ct_max_value[i][j].size(); ++k) {
                if (gt(ct_max_value[i][j][k], out_max_value[i][j])) {
                    out_max_value[i][j] = ct_max_value[i][j][k];
                }
            }
        }
    }
}

template <typename DCRTPoly>
void LibsnarkProofSystem::ConstrainKeySwitch(const Ciphertext<DCRTPoly>& ctxt_in, const EvalKey<DCRTPoly>& ek,
                                             Ciphertext<DCRTPoly>& ctxt_out) {
    assert(ctxt_in->GetElements().size() == 2);
    auto in = *GetProofMetadata(ctxt_in);

    std::vector<DCRTPoly>& cv = ctxt_in->GetElements();

    //    std::shared_ptr<std::vector<DCRTPoly>> ba = KeySwitchCore(cv[1], ek) ;
    auto a                                        = cv[1];
    const auto cryptoParamsBase                   = ek->GetCryptoParameters();
    std::shared_ptr<std::vector<DCRTPoly>> digits = KeySwitchBV().EvalKeySwitchPrecomputeCore(a, cryptoParamsBase);
    vector<vector<vector<pb_linear_combination<FieldT>>>> digits_lc;
    vector<vector<FieldT>> digits_max_value;
    ConstrainKeySwitchPrecomputeCore(a, cryptoParamsBase, digits, in[1], in.max_value[1], digits_lc, digits_max_value);

    std::shared_ptr<std::vector<DCRTPoly>> result = KeySwitchBV().EvalFastKeySwitchCore(digits, ek, a.GetParams());
    vector<vector<vector<pb_linear_combination<FieldT>>>> out_lc;
    vector<vector<FieldT>> out_max_value;
    ConstrainFastKeySwitchCore(digits, ek, a.GetParams(), result, digits_lc, digits_max_value, out_lc, out_max_value);
    LibsnarkProofMetadata out(out_lc);
    out.max_value = out_max_value;
    out.modulus   = out.modulus;

    auto ba = digits;
    assert(cv[0].GetFormat() == (*ba)[0].GetFormat());
    //    cv[0] += (*ba)[0];
    for (size_t j = 0; j < out[0].size(); j++) {
        for (size_t k = 0; k < out[0][j].size(); k++) {
            LazyAddModGadget<FieldT> g(pb, out[0][j][k], out.max_value[0][j], in[0][j][k], in.max_value[0][j],
                                       in.modulus[j]);
        }
    }

    assert(cv[1].GetFormat() == (*ba)[1].GetFormat());
    //    if (cv.size() > 2) {
    //        cv[1] += (*ba)[1];
    //    }
    //    else {
    cv[1] = (*ba)[1];
    //    }

    SetProofMetadata(ctxt_out, std::make_shared<LibsnarkProofMetadata>(out));
}

template <typename DCRTPoly>
void LibsnarkProofSystem::ConstrainRelin(const Ciphertext<DCRTPoly>& ciphertext, Ciphertext<DCRTPoly>& out) {
    assert(!!ciphertext);
    assert(!!out);

    const LibsnarkProofMetadata in = *(GetProofMetadata(ciphertext));
    const size_t num_poly          = ciphertext->GetElements().size();
    const size_t num_limbs         = ciphertext->GetElements()[0].GetNumOfElements();
    assert(num_poly == 3);  // We don't support higher-order relin
    assert(in.size() == num_poly);
    assert(in[0].size() == num_limbs);
    assert(out->GetElements()[0].GetNumOfElements() == num_limbs);

    LibsnarkProofMetadata out_metadata(2);

    const auto evalKeyVec = ciphertext->GetCryptoContext()->GetEvalMultKeyVector(ciphertext->GetKeyTag());
    assert(evalKeyVec.size() >= (ciphertext->GetElements().size() - 2));

    const auto cv = ciphertext->GetElements();
    for (auto& c : cv) {
        assert(c.GetFormat() == Format::EVALUATION);  // Should always hold for BGV, no need to constrain
        // c.SetFormat(Format::EVALUATION);
    }

    auto algo = ciphertext->GetCryptoContext()->GetScheme();

    out_metadata           = LibsnarkProofMetadata(2);
    out_metadata[0]        = in[0];
    out_metadata[1]        = in[1];
    out_metadata.max_value = {in.max_value[0], in.max_value[1]};
    out_metadata.modulus   = in.modulus;

    for (size_t j = 2; j < num_poly; j++) {
        //        auto ab      = algo->KeySwitchCore(cv[j], evalKeyVec[j - 2]);
        auto evalKey = evalKeyVec[j - 2];

        const auto cryptoParamsBase = evalKey->GetCryptoParameters();
        std::shared_ptr<std::vector<DCRTPoly>> digits =
            KeySwitchBV().EvalKeySwitchPrecomputeCore(cv[j], cryptoParamsBase);

        vector<vector<vector<pb_linear_combination<FieldT>>>> tmp_lc;
        vector<vector<FieldT>> tmp_max_value;
        ConstrainKeySwitchPrecomputeCore(cv[j], evalKey->GetCryptoParameters(), digits, in[j], in.max_value[j], tmp_lc,
                                         tmp_max_value);

        std::shared_ptr<std::vector<DCRTPoly>> result =
            KeySwitchBV().EvalFastKeySwitchCore(digits, evalKey, cv[j].GetParams());
        vector<vector<vector<pb_linear_combination<FieldT>>>> tmp2_lc;
        vector<vector<FieldT>> tmp2_max_value;
        ConstrainFastKeySwitchCore(digits, evalKey, cv[j].GetParams(), result, tmp_lc, tmp_max_value, tmp2_lc,
                                   tmp2_max_value);

        //        cv[0] += (*ab)[0];
        //        cv[1] += (*ab)[1];
        LibsnarkProofMetadata tmp_metadata(tmp2_lc);
        tmp_metadata.max_value = tmp2_max_value;
        tmp_metadata.modulus   = in.modulus;
        constrain_addmod_lazy(out_metadata, 0, tmp_metadata, 0, out_metadata, 0);
        constrain_addmod_lazy(out_metadata, 1, tmp_metadata, 1, out_metadata, 1);
    }

    SetProofMetadata(out, std::make_shared<LibsnarkProofMetadata>(out_metadata));
}

template <typename DCRTPoly>
void LibsnarkProofSystem::ConstrainRotate(const Ciphertext<DCRTPoly>& ciphertext, const int rot_idx,
                                          Ciphertext<DCRTPoly>& ctxt_out) {
    auto in = *GetProofMetadata(ciphertext);

    auto index    = rot_idx;
    const auto cc = ciphertext->GetCryptoContext();

    auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(ciphertext->GetKeyTag());
    //    return GetScheme()->EvalAtIndex(ciphertext, index, evalKeyMap);

    usint M = ciphertext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder();

    auto algo          = ciphertext->GetCryptoContext()->GetScheme();
    uint32_t autoIndex = FindAutomorphismIndex2n(index, M);  // TODO: use scheme stored in ctxt instead

    //    return EvalAutomorphism(ciphertext, autoIndex, evalKeyMap);

    auto i   = autoIndex;
    auto key = evalKeyMap.find(i);
    assert(key != evalKeyMap.end());

    //    if (key == evalKeyMap.end()) {
    //        std::string errorMsg(std::string("Could not find an EvalKey for index ") + std::to_string(i) + CALLER_INFO);
    //        OPENFHE_THROW(type_error, errorMsg);
    //    }

    auto evalKey = key->second;

    //    CheckKey(evalKey);

    //    return GetScheme()->EvalAutomorphism(ciphertext, i, evalKeyMap);

    auto evalKeyIterator = evalKeyMap.find(i);
    assert(evalKeyIterator != evalKeyMap.end());
    using Element                  = DCRTPoly;
    const std::vector<Element>& cv = ciphertext->GetElements();
    usint N                        = cv[0].GetRingDimension();

    std::vector<usint> vec(N);
    PrecomputeAutoMap(N, i, &vec);

    //    auto algo = ciphertext->GetCryptoContext()->GetScheme();

    Ciphertext<Element> result = ciphertext->Clone();

    auto old(result);
    algo->KeySwitchInPlace(result, evalKeyIterator->second);
    ConstrainKeySwitch(old, evalKeyIterator->second, result);

    //    rcv[0] = rcv[0].AutomorphismTransform(i, vec);
    //    rcv[1] = rcv[1].AutomorphismTransform(i, vec);
    auto precomp = vec;
    LibsnarkProofMetadata out(in);
    for (size_t i = 0; i < out.size(); i++) {
        for (size_t j = 0; j < out[0].size(); ++j) {
            for (size_t k = 0; k < out[0][0].size(); ++k) {
                out[i][j][k] = out[i][j][precomp[k]];
            }
        }
    }

    SetProofMetadata(ctxt_out, std::make_shared<LibsnarkProofMetadata>(out));
}

#endif  //OPENFHE_PROOFSYSTEM_LIBSNARK_CPP
