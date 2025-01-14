#ifndef OPENFHE_PROOFSYSTEM_LIBSNARK_H
#define OPENFHE_PROOFSYSTEM_LIBSNARK_H

#include "proofsystem.h"
#include "libff/algebra/field_utils/field_utils.hpp"
#include "libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp"
#include "libsnark/common/default_types/r1cs_ppzksnark_pp.hpp"
#include "libsnark/gadgetlib1/pb_variable.hpp"
#include "gadgets_libsnark.h"
#include "key/evalkey.h"
#include <libff/algebra/fields/binary/gf192.hpp>

#include <vector>
using std::vector;

using namespace libsnark;
// typedef libff::Fr<default_r1cs_ppzksnark_pp> FieldT;
typedef libff::gf192 FieldT;

class LibsnarkWitnessMetadata : public ProofsystemMetadata {
public:
    LibsnarkWitnessMetadata() = default;
    size_t index_start;
    std::shared_ptr<protoboard<FieldT>> pb;
    vector<std::shared_ptr<gadget_gen<FieldT>>> gadgets;

    static std::string GetKey() {
        return std::string("proofsystem_metadata_") + typeid(LibsnarkWitnessMetadata).name();
    };

    template <typename GadgetT>
    std::enable_if_t<std::is_base_of_v<gadget_gen<FieldT>, GadgetT>, void> add_gadget(const GadgetT& gadget) {
        gadgets.push_back(std::make_shared<GadgetT>(gadget));
    }

    template <typename GadgetT>
    std::enable_if_t<std::is_base_of_v<gadget_gen<FieldT>, GadgetT>, void> add_gadget(GadgetT&& gadget) {
        gadgets.push_back(std::make_shared<GadgetT>(gadget));
    }
};

class LibsnarkConstraintMetadata : public ProofsystemMetadata,
                                   private vector<vector<vector<pb_linear_combination<FieldT>>>> {
public:
    WireID wire_id;
    LibsnarkWitnessMetadata witness_metadata;
    vector<size_t> modulus;
    vector<vector<FieldT>> max_value;

    explicit LibsnarkConstraintMetadata(size_t n = 0)
        : ProofsystemMetadata(), vector<vector<vector<pb_linear_combination<FieldT>>>>(n), modulus(0), max_value(n) {}

    explicit LibsnarkConstraintMetadata(
        const vector<vector<vector<pb_linear_combination<FieldT>>>>& pb_linear_combinations)
        : ProofsystemMetadata(),
          vector<vector<vector<pb_linear_combination<FieldT>>>>(pb_linear_combinations),
          modulus(pb_linear_combinations.size()),
          max_value(pb_linear_combinations.size()) {}

    using vector<vector<vector<pb_linear_combination<FieldT>>>>::vector;
    using vector<vector<vector<pb_linear_combination<FieldT>>>>::operator[];
    using vector<vector<vector<pb_linear_combination<FieldT>>>>::at;
    using vector<vector<vector<pb_linear_combination<FieldT>>>>::size;
    using vector<vector<vector<pb_linear_combination<FieldT>>>>::operator=;
    using vector<vector<vector<pb_linear_combination<FieldT>>>>::push_back;
    using vector<vector<vector<pb_linear_combination<FieldT>>>>::emplace_back;

    static std::string GetKey() {
        return std::string("proofsystem_metadata_") + typeid(LibsnarkConstraintMetadata).name();
    };

    inline size_t get_bit_size(size_t i, size_t j) const {
        return max_value[i][j].as_bigint().num_bits();
    }

    inline std::tuple<size_t, size_t, size_t> get_dims() const {
        return std::make_tuple(size(), this->operator[](0).size(), this->operator[](0)[0].size());
    }

    inline bool matches_dim(const LibsnarkConstraintMetadata& other) const {
        auto [n1, n2, n3] = get_dims();
        auto [m1, m2, m3] = other.get_dims();
        return (n1 == m1) && (n2 == m2) && (n3 == m3);
    }

    template <typename Element>
    inline bool matches_dim(ConstCiphertext<Element>& ciphertext) const {
        auto [n1, n2, n3] = get_dims();
        auto [m1, m2, m3] =
            std::make_tuple(ciphertext->GetElements().size(), ciphertext->GetElements()[0].GetNumOfElements(),
                            ciphertext->GetElements()[0].GetElementAtIndex(0).GetLength());
        return (n1 == m1) && (n2 == m2) && (n3 == m3);
    }
};

class LibsnarkProofSystem : public ProofSystem<DCRTPoly, LibsnarkConstraintMetadata, LibsnarkWitnessMetadata> {
protected:
    vector<std::shared_ptr<gadget_gen<FieldT>>> constrain_addmod_lazy(const LibsnarkConstraintMetadata& in1,
                                                                      size_t index_1,
                                                                      const LibsnarkConstraintMetadata& in2,
                                                                      size_t index_2, LibsnarkConstraintMetadata& out,
                                                                      size_t index_out);
    void constrain_addmod_lazy(const LibsnarkConstraintMetadata& in1, size_t index_1,
                               const LibsnarkConstraintMetadata& in2, size_t index_2, LibsnarkConstraintMetadata& out,
                               size_t index_out, vector<std::shared_ptr<gadget_gen<FieldT>>>& gadgets_append);
    void constrain_submod_lazy(const LibsnarkConstraintMetadata& in1, size_t index_1,
                               const LibsnarkConstraintMetadata& in2, size_t index_2, LibsnarkConstraintMetadata& out,
                               size_t index_out);
    vector<std::shared_ptr<gadget_gen<FieldT>>> constrain_mulmod_lazy(const LibsnarkConstraintMetadata& in1,
                                                                      size_t index_1,
                                                                      const LibsnarkConstraintMetadata& in2,
                                                                      size_t index_2, LibsnarkConstraintMetadata& out,
                                                                      size_t index_out);
    void constrain_mulmod_lazy(const LibsnarkConstraintMetadata& in1, size_t index_1,
                               const LibsnarkConstraintMetadata& in2, size_t index_2, LibsnarkConstraintMetadata& out,
                               size_t index_out, vector<std::shared_ptr<gadget_gen<FieldT>>>& gadgets_append);
    std::unordered_map<WireID, LibsnarkWitnessMetadata> wire_metadata;
    size_t GetGlobalWireId() override;
    void SetGlobalWireId(size_t globalWireId) override;

    size_t global_wire_id = 0;

public:
    protoboard<FieldT> pb;

    explicit LibsnarkProofSystem(const CryptoContext<DCRTPoly>& cc)
        : ProofSystem<DCRTPoly, LibsnarkConstraintMetadata, LibsnarkWitnessMetadata>(cc) {
        default_r1cs_ppzksnark_pp::init_public_params();
    }

    void SetMode(PROOFSYSTEM_MODE mode) override {
        ProofSystem<DCRTPoly, LibsnarkConstraintMetadata, LibsnarkWitnessMetadata>::SetMode(mode);
        SetGlobalWireId(0);
        if (mode == PROOFSYSTEM_MODE_CONSTRAINT_GENERATION) {
            pb = protoboard<FieldT>();
        }
    }
    size_t GetNextWireId() override;

    std::shared_ptr<LibsnarkConstraintMetadata> ConstrainPublicOutput(Ciphertext<DCRTPoly>& ciphertext);

    LibsnarkConstraintMetadata PublicInputConstraint(ConstCiphertext<DCRTPoly> ciphertext) override;
    void PublicInputWitness(ConstCiphertext<DCRTPoly> ciphertext) override;

    LibsnarkConstraintMetadata EvalAddConstraint(const LibsnarkConstraintMetadata& m1,
                                                 const LibsnarkConstraintMetadata& m2) override;
    void EvalAddWitness(ConstCiphertext<DCRTPoly> ciphertext1, ConstCiphertext<DCRTPoly> ciphertext2,
                        ConstCiphertext<DCRTPoly> ciphertext_out) override;

    LibsnarkConstraintMetadata EvalSubConstraint(const LibsnarkConstraintMetadata& m1,
                                                 const LibsnarkConstraintMetadata& m2) override;
    void EvalSubWitness(ConstCiphertext<DCRTPoly> ciphertext1, ConstCiphertext<DCRTPoly> ciphertext2,
                        ConstCiphertext<DCRTPoly> ciphertext_out) override;

    LibsnarkConstraintMetadata EvalMultNoRelinConstraint(const LibsnarkConstraintMetadata& m1,
                                                         const LibsnarkConstraintMetadata& m2) override;
    void EvalMultNoRelinWitness(ConstCiphertext<DCRTPoly> ciphertext1, ConstCiphertext<DCRTPoly> ciphertext2,
                                ConstCiphertext<DCRTPoly> ciphertext_out) override;

    LibsnarkConstraintMetadata EvalSquareConstraint(const LibsnarkConstraintMetadata& m) override;
    void EvalSquareWitness(ConstCiphertext<DCRTPoly> ciphertext, ConstCiphertext<DCRTPoly> ciphertext_out) override;

    LibsnarkConstraintMetadata RescaleConstraint(ConstCiphertext<DCRTPoly> ctxt_in,
                                                 const LibsnarkConstraintMetadata& m) override;
    void RescaleWitness(ConstCiphertext<DCRTPoly> ciphertext, ConstCiphertext<DCRTPoly> ciphertext_out) override;

    LibsnarkConstraintMetadata EvalRotateConstraint(ConstCiphertext<DCRTPoly> ciphertext, int rot_idx,
                                                    ConstCiphertext<DCRTPoly> ctxt_out) override;
    void EvalRotateWitness(ConstCiphertext<DCRTPoly> ciphertext, int rot_idx,
                           ConstCiphertext<DCRTPoly> ctxt_out) override;

    LibsnarkConstraintMetadata RelinearizeConstraint(ConstCiphertext<DCRTPoly> ctxt_in) override;
    void RelinearizeWitness(ConstCiphertext<DCRTPoly> ciphertext, ConstCiphertext<DCRTPoly> ciphertext_out) override;

    LibsnarkConstraintMetadata KeySwitchConstraint(ConstCiphertext<DCRTPoly> ctxt_in, const EvalKey<DCRTPoly>& ek,
                                                   ConstCiphertext<DCRTPoly> ctxt_out) override;
    void KeySwitchWitness(ConstCiphertext<DCRTPoly> ctxt_in, const EvalKey<DCRTPoly>& ek,
                          ConstCiphertext<DCRTPoly>& ctxt_out) override;

    void EncryptWitness(Plaintext plaintext, PublicKey<DCRTPoly> publicKey,
                        ConstCiphertext<DCRTPoly> ciphertext) override;
    void EncryptConstraint(Plaintext plaintext, DoublePublicKey<DCRTPoly> publicKey,
                           DoubleCiphertext<DCRTPoly> ciphertext) override;
    void EncryptWitness(Plaintext plaintext, DoublePublicKey<DCRTPoly> publicKey,
                        DoubleCiphertext<DCRTPoly> ciphertext) override;

    LibsnarkConstraintMetadata EncryptConstraint(Plaintext plaintext, PublicKey<DCRTPoly> publicKey) override;

    void ConstrainSquare2(ConstCiphertext<DCRTPoly>& ctxt, Ciphertext<DCRTPoly>& ctxt_out);

    void ConstrainRescale(ConstCiphertext<DCRTPoly>& ctxt, Ciphertext<DCRTPoly>& ctxt_out);

    void SetFormatConstraint(const Format format, const DCRTPoly::PolyType& in, const DCRTPoly::PolyType& out,
                             const vector<pb_linear_combination<FieldT>>& in_lc, const FieldT& in_max_value,
                             vector<pb_linear_combination<FieldT>>& out_lc, FieldT& out_max_value,
                             LibsnarkWitnessMetadata& witness_metadata);

    void SetFormatConstraint(const Format format, const DCRTPoly& in, const DCRTPoly& out,
                             const vector<vector<pb_linear_combination<FieldT>>>& in_lc,
                             const vector<FieldT>& in_max_value, vector<vector<pb_linear_combination<FieldT>>>& out_lc,
                             vector<FieldT>& out_max_value, LibsnarkWitnessMetadata& witness_metadata);

    void SetFormatWitness(const Format format, LibsnarkWitnessMetadata& witness_metadata);

    void NTTOpenfheConstraint(const DCRTPoly::PolyType::Vector& rootOfUnityTable,
                              const DCRTPoly::PolyType::Vector& preconRootOfUnityTable,
                              const DCRTPoly::PolyType& element_in, const DCRTPoly::PolyType& element_out,
                              const vector<pb_linear_combination<FieldT>>& in_lc, const FieldT& in_max_value,
                              vector<pb_linear_combination<FieldT>>& out_lc, FieldT& out_max_value,
                              LibsnarkWitnessMetadata& witness_metadata);

    void NTTOpenfheWitness(LibsnarkWitnessMetadata& witness_metadata);

    void NTTLinalgConstraint(const DCRTPoly::PolyType& element_in, const DCRTPoly::PolyType& element_out,
                             const vector<pb_linear_combination<FieldT>>& in_lc, const FieldT& in_max_value,
                             vector<pb_linear_combination<FieldT>>& out_lc, FieldT& out_max_value,
                             LibsnarkWitnessMetadata& witness_metadata);

    void NTTLinalgWitness(LibsnarkWitnessMetadata& witness_metadata);

    void INTTOpenfheConstraint(const DCRTPoly::PolyType::Vector& rootOfUnityInverseTable,
                               const DCRTPoly::PolyType::Vector& preconRootOfUnityInverseTable,
                               const DCRTPoly::PolyType::Vector::Integer& cycloOrderInv,
                               const DCRTPoly::PolyType::Vector::Integer& preconCycloOrderInv,
                               const DCRTPoly::PolyType& element_in, const DCRTPoly::PolyType& element_out,
                               const vector<pb_linear_combination<FieldT>>& in_lc, const FieldT& in_max_value,
                               vector<pb_linear_combination<FieldT>>& out_lc, FieldT& out_max_value,
                               LibsnarkWitnessMetadata& witness_metadata);

    void INTTOpenfheWitness(LibsnarkWitnessMetadata& witness_metadata);

    void SwitchModulusConstraint(const DCRTPoly::PolyType::Vector::Integer& newModulus,
                                 const DCRTPoly::PolyType::Vector::Integer& rootOfUnity, const DCRTPoly::PolyType& in,
                                 const DCRTPoly::PolyType& out, const vector<pb_linear_combination<FieldT>>& in_lc,
                                 const FieldT& in_max_value, vector<pb_linear_combination<FieldT>>& out_lc,
                                 FieldT& out_max_value, LibsnarkWitnessMetadata& witness_metadata);

    void SwitchModulusWitness(LibsnarkWitnessMetadata& witness_metadata);

    void ConstrainKeySwitchPrecomputeCore(
        const DCRTPoly& in, const std::shared_ptr<CryptoParametersBase<DCRTPoly>>& cryptoParamsBase,
        const std::shared_ptr<std::vector<DCRTPoly>>& out, const vector<vector<pb_linear_combination<FieldT>>>& in_lc,
        const vector<FieldT>& in_max_value, vector<vector<vector<pb_linear_combination<FieldT>>>>& out_lc,
        vector<vector<FieldT>>& out_max_value, LibsnarkWitnessMetadata& witness_metadata);

    void ConstrainFastKeySwitchCore(const EvalKey<DCRTPoly>& evalKey, const std::shared_ptr<DCRTPoly::Params>& paramsQl,
                                    const vector<vector<vector<pb_linear_combination<FieldT>>>>& in_lc,
                                    const vector<vector<FieldT>>& in_max_value,
                                    vector<vector<vector<pb_linear_combination<FieldT>>>>& out_lc,
                                    vector<vector<FieldT>>& out_max_value);

    void ConstrainFastKeySwitchCore(const std::shared_ptr<std::vector<DCRTPoly>>& digits,
                                    const EvalKey<DCRTPoly>& evalKey, const std::shared_ptr<DCRTPoly::Params>& paramsQl,
                                    std::shared_ptr<std::vector<DCRTPoly>>& out,
                                    const vector<vector<vector<pb_linear_combination<FieldT>>>>& in_lc,
                                    const vector<vector<FieldT>>& in_max_value,
                                    vector<vector<vector<pb_linear_combination<FieldT>>>>& out_lc,
                                    vector<vector<FieldT>>& out_max_value);

    void FinalizeOutputConstraints(Ciphertext<DCRTPoly>& ctxt, const LibsnarkConstraintMetadata& out_vars) override;
};

void NTTParameters(DCRTPoly::PolyType::Vector::Integer rootOfUnity, usint CycloOrder,
                   DCRTPoly::PolyType::Integer modulus, DCRTPoly::PolyType::Vector& rootOfUnityTable,
                   DCRTPoly::PolyType::Vector& preconRootOfUnityTable);

void INTTParameters(DCRTPoly::PolyType::Vector::Integer rootOfUnity, usint CycloOrder,
                    DCRTPoly::PolyType::Integer modulus, DCRTPoly::PolyType::Vector& rootOfUnityInverseTable,
                    DCRTPoly::PolyType::Vector& preconRootOfUnityInverseTable,
                    DCRTPoly::PolyType::Vector::Integer& cycloOrderInv,
                    DCRTPoly::PolyType::Vector::Integer& preconCycloOrderInv);

#endif  //OPENFHE_PROOFSYSTEM_LIBSNARK_H
