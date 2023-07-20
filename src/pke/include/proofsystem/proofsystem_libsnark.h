#ifndef OPENFHE_PROOFSYSTEM_LIBSNARK_H
#define OPENFHE_PROOFSYSTEM_LIBSNARK_H

#include "proofsystem.h"
#include "libff/algebra/fields/field_utils.hpp"
#include "libsnark/zk_proof_systems/ppzksnark/r1cs_ppzksnark/r1cs_ppzksnark.hpp"
#include "libsnark/common/default_types/r1cs_ppzksnark_pp.hpp"
#include "libsnark/gadgetlib1/pb_variable.hpp"
#include <vector>
using std::vector;

using namespace libsnark;
typedef libff::Fr<default_r1cs_ppzksnark_pp> FieldT;

class LibsnarkProofMetadata : public ProofMetadata, private vector<vector<vector<pb_linear_combination<FieldT>>>> {
public:
    explicit LibsnarkProofMetadata(const vector<vector<vector<pb_linear_combination<FieldT>>>>& pb_linear_combinations)
        : ProofMetadata(), vector<vector<vector<pb_linear_combination<FieldT>>>>(pb_linear_combinations) {}

    using vector<vector<vector<pb_linear_combination<FieldT>>>>::vector;
    using vector<vector<vector<pb_linear_combination<FieldT>>>>::operator[];
    using vector<vector<vector<pb_linear_combination<FieldT>>>>::size;
    using vector<vector<vector<pb_linear_combination<FieldT>>>>::operator=;
};

class LibsnarkProofSystem : ProofSystem {
public:
    protoboard<FieldT> pb;

    LibsnarkProofSystem() {
        default_r1cs_ppzksnark_pp::init_public_params();
        pb = protoboard<FieldT>();
        // TODO FIXME this is a hack
        pb_variable<FieldT> dummy;
        dummy.allocate(pb);
        pb.set_input_sizes(1);
    }
    void ConstrainPublicInput(Ciphertext<DCRTPoly>& ciphertext) override;
    void ConstrainAddition(const Ciphertext<DCRTPoly>& ctxt1, const Ciphertext<DCRTPoly>& ctxt2,
                           Ciphertext<DCRTPoly>& ctxt_out) override;
    static std::shared_ptr<LibsnarkProofMetadata> GetProofMetadata(const Ciphertext<DCRTPoly>& ciphertext);
    static void SetProofMetadata(const Ciphertext<DCRTPoly>& ciphertext, const std::shared_ptr<LibsnarkProofMetadata>& metadata);

};
#endif  //OPENFHE_PROOFSYSTEM_LIBSNARK_H