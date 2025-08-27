#pragma once
#ifndef AAGOS_MUT_LANDSCAPE_INFO_H
#define AAGOS_MUT_LANDSCAPE_INFO_H

#include "AagosOrg.hpp"

#include "emp/Evolve/Systematics.hpp"

namespace aagos {

struct AagosMutLandscapeInfo : public emp::datastruct::mut_landscape_info<AagosOrg::Phenotype> {

  bool HasMutationType(const std::string& mut_type) const {
    return emp::Has(mut_counts, mut_type);
  }

  int GetMutationCount(const std::string& mut_type) const {
    emp_assert(HasMutationType(mut_type));
    return mut_counts.at(mut_type);
  }

};

}

#endif