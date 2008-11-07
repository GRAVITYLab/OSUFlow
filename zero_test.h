#define fuzz_factor 1.0e-7

#define Zero_test1(result,term1) \
  if (fabs(result) <= fuzz_factor * fabs(term1)) \
    result = 0.0

#define Zero_test2(result,term1,term2) \
  if (fabs(result) <= fuzz_factor * (fabs(term1) + fabs(term2))) \
    result = 0.0

#define Zero_test3(result,term1,term2,term3) \
  if (fabs(result) <= fuzz_factor * (fabs(term1) + fabs(term2) + fabs(term3))) \
    result = 0.0

#define Zero_test4(result,term1,term2,term3, term4) \
  if (fabs(result) <= fuzz_factor * (fabs(term1) + fabs(term2) + fabs(term3) + fabs(term4))) \
    result = 0.0

#define Zero_test6(result,term1,term2,term3, term4, term5, term6) \
  if (fabs(result) <= fuzz_factor * (fabs(term1) + fabs(term2) + fabs(term3) + fabs(term4) + fabs(term5) + fabs(term6))) \
    result = 0.0
