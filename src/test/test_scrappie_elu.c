#include <CUnit/Basic.h>
#include <stdbool.h>

#include <util.h>

/**  Initialise test
 *
 *   @returns 0 on success, non-zero on failure
 **/
int init_test_elu(void) {
    return 0;
}

/**  Clean up after test
 *
 *   @returns 0 on success, non-zero on failure
 **/
int clean_test_elu(void) {
    return 0;
}

void test_zero_elu(void) {
    __m128 x = elufv(_mm_setzero_ps());
    __v4sf cmp = _mm_cmpneq_ps(x, _mm_setzero_ps());

    CU_ASSERT_EQUAL(cmp[0], 0.0f);
    CU_ASSERT_EQUAL(cmp[1], 0.0f);
    CU_ASSERT_EQUAL(cmp[2], 0.0f);
    CU_ASSERT_EQUAL(cmp[3], 0.0f);
}

void test_negzero_elu(void) {
    __m128 x = elufv(-_mm_setzero_ps());
    __v4sf cmp = _mm_cmpneq_ps(x, _mm_setzero_ps());

    CU_ASSERT_EQUAL(cmp[0], 0.0f);
    CU_ASSERT_EQUAL(cmp[1], 0.0f);
    CU_ASSERT_EQUAL(cmp[2], 0.0f);
    CU_ASSERT_EQUAL(cmp[3], 0.0f);
}

void test_positive_elu(void) {
    __m128 x = _mm_setr_ps(1.0, 2.0, 3.0, 4.0);
    __v4sf y = elufv(x);

    CU_ASSERT_EQUAL(y[0], 1.0f);
    CU_ASSERT_EQUAL(y[1], 2.0f);
    CU_ASSERT_EQUAL(y[2], 3.0f);
    CU_ASSERT_EQUAL(y[3], 4.0f);
}

void test_negative_elu(void) {
    __m128 x = _mm_setr_ps(-1.0, -2.0, -3.0, -4.0);
    __v4sf y = elufv(x);

    CU_ASSERT_DOUBLE_EQUAL(y[0], -0.6321206f, 1e-6f);
    CU_ASSERT_DOUBLE_EQUAL(y[1], -0.8646647f, 1e-6f);
    CU_ASSERT_DOUBLE_EQUAL(y[2], -0.9502129f, 1e-6f);
    CU_ASSERT_DOUBLE_EQUAL(y[3], -0.9816844f, 1e-6f);
}

void test_mixed_elu(void) {
    __m128 x = _mm_setr_ps(1.0, -2.0, 3.0, -4.0);
    __v4sf y = elufv(x);

    CU_ASSERT_EQUAL(y[0], 1.0f);
    CU_ASSERT_DOUBLE_EQUAL(y[1], -0.8646647f, 1e-6f);
    CU_ASSERT_EQUAL(y[2], 3.0f);
    CU_ASSERT_DOUBLE_EQUAL(y[3], -0.9816844f, 1e-6f);
}

/**   Register tests with CUnit
 *
 *    @returns 0 on success, non-zero on failure
 **/
int register_test_elu(void) {
    CU_pSuite suite = CU_add_suite("Test output of vectorised ELU unit",
                                   init_test_elu,
                                   clean_test_elu);
    if (NULL == suite) {
        return CU_get_error();
    }

    if (NULL ==
        CU_add_test(suite, "Test input all zeros", test_zero_elu)) {
        return CU_get_error();
    }
    if (NULL ==
        CU_add_test(suite, "Test input all negative zeros", test_negzero_elu)) {
        return CU_get_error();
    }
    if (NULL ==
        CU_add_test(suite, "Test input all positive", test_positive_elu)) {
        return CU_get_error();
    }
    if (NULL ==
        CU_add_test(suite, "Test input all negative", test_negative_elu)) {
        return CU_get_error();
    }
    if (NULL ==
        CU_add_test(suite, "Test input mixed positive and negative", test_mixed_elu)) {
        return CU_get_error();
    }

    return 0;
}
