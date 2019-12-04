set(FF_EXTRASRCS)
set(FF_EXTRALIBS)
if(${CURVE} STREQUAL "BN128")
  set(
    FF_EXTRASRCS

    ${FF_EXTRASRCS}
    algebra/curves/bn128/bn128_g1.cpp
    algebra/curves/bn128/bn128_g2.cpp
    algebra/curves/bn128/bn128_gt.cpp
    algebra/curves/bn128/bn128_init.cpp
    algebra/curves/bn128/bn128_pairing.cpp
    algebra/curves/bn128/bn128_pp.cpp
  )
  set(
    FF_EXTRALIBS

    ${FF_EXTRALIBS}
    zm
  )
endif()

add_library(
  ff
  STATIC

  common/double.cpp
  common/profiling.cpp
  common/utils.cpp

  algebra/curves/toy_curve/toy_curve_g1.cpp
  algebra/curves/toy_curve/toy_curve_g2.cpp
  algebra/curves/toy_curve/toy_curve_init.cpp
  algebra/curves/toy_curve/toy_curve_pairing.cpp
  algebra/curves/toy_curve/toy_curve_pp.cpp

  algebra/curves/alt_bn128/alt_bn128_g1.cpp
  algebra/curves/alt_bn128/alt_bn128_g2.cpp
  algebra/curves/alt_bn128/alt_bn128_init.cpp
  algebra/curves/alt_bn128/alt_bn128_pairing.cpp
  algebra/curves/alt_bn128/alt_bn128_pp.cpp

  algebra/curves/bw12_446/bw12_446_g1.cpp
  algebra/curves/bw12_446/bw12_446_g2.cpp
  algebra/curves/bw12_446/bw12_446_init.cpp
  algebra/curves/bw12_446/bw12_446_pairing.cpp
  algebra/curves/bw12_446/bw12_446_pp.cpp

  algebra/curves/edwards/edwards_g1.cpp
  algebra/curves/edwards/edwards_g2.cpp
  algebra/curves/edwards/edwards_init.cpp
  algebra/curves/edwards/edwards_pairing.cpp
  algebra/curves/edwards/edwards_pp.cpp

  algebra/curves/mnt/mnt46_common.cpp

  algebra/curves/mnt/mnt4/mnt4_g1.cpp
  algebra/curves/mnt/mnt4/mnt4_g2.cpp
  algebra/curves/mnt/mnt4/mnt4_init.cpp
  algebra/curves/mnt/mnt4/mnt4_pairing.cpp
  algebra/curves/mnt/mnt4/mnt4_pp.cpp

  algebra/curves/mnt/mnt6/mnt6_g1.cpp
  algebra/curves/mnt/mnt6/mnt6_g2.cpp
  algebra/curves/mnt/mnt6/mnt6_init.cpp
  algebra/curves/mnt/mnt6/mnt6_pairing.cpp
  algebra/curves/mnt/mnt6/mnt6_pp.cpp

  algebra/curves/mnt753/mnt46753_common.cpp

  algebra/curves/mnt753/mnt6753/mnt6753_g1.cpp
  algebra/curves/mnt753/mnt6753/mnt6753_g2.cpp
  algebra/curves/mnt753/mnt6753/mnt6753_init.cpp
  algebra/curves/mnt753/mnt6753/mnt6753_pairing.cpp
  algebra/curves/mnt753/mnt6753/mnt6753_pp.cpp

  algebra/curves/mnt753/mnt4753/mnt4753_g1.cpp
  algebra/curves/mnt753/mnt4753/mnt4753_g2.cpp
  algebra/curves/mnt753/mnt4753/mnt4753_init.cpp
  algebra/curves/mnt753/mnt4753/mnt4753_pairing.cpp
  algebra/curves/mnt753/mnt4753/mnt4753_pp.cpp

  algebra/curves/pendulum/pendulum_g1.cpp
  algebra/curves/pendulum/pendulum_g2.cpp
  algebra/curves/pendulum/pendulum_init.cpp
  algebra/curves/pendulum/pendulum_pairing.cpp
  algebra/curves/pendulum/pendulum_pp.cpp

  algebra/curves/bls12_377/bls12_377_g1.cpp
  algebra/curves/bls12_377/bls12_377_g2.cpp
  algebra/curves/bls12_377/bls12_377_init.cpp
  algebra/curves/bls12_377/bls12_377_pairing.cpp
  algebra/curves/bls12_377/bls12_377_pp.cpp

  algebra/curves/sw6/sw6_g1.cpp
  algebra/curves/sw6/sw6_g2.cpp
  algebra/curves/sw6/sw6_init.cpp
  algebra/curves/sw6/sw6_pairing.cpp
  algebra/curves/sw6/sw6_pp.cpp

  ${FF_EXTRASRCS}
)
target_link_libraries(
  ff

  GMP::gmp
  ${PROCPS_LIBRARIES}
  ${FF_EXTRALIBS}
)
target_include_directories(
  ff
  PUBLIC ..
)

install(
  DIRECTORY "" DESTINATION "include/libff"
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.tcc"
  PATTERN "tests" EXCLUDE
  PATTERN "examples" EXCLUDE
)

install(
  TARGETS ff DESTINATION lib
)

# Tests
if ("${IS_LIBFF_PARENT}")
  add_executable(
    algebra_bilinearity_test
    EXCLUDE_FROM_ALL

    algebra/curves/tests/test_bilinearity.cpp
  )
  target_link_libraries(
    algebra_bilinearity_test

    ff
  )

  add_executable(
    algebra_groups_test
    EXCLUDE_FROM_ALL

    algebra/curves/tests/test_groups.cpp
  )
  target_link_libraries(
    algebra_groups_test

    ff
  )

  add_executable(
    algebra_fields_test
    EXCLUDE_FROM_ALL

    algebra/fields/tests/test_fields.cpp
  )
  target_link_libraries(
    algebra_fields_test

    ff
  )

  include(CTest)
  add_test(
    NAME algebra_bilinearity_test
    COMMAND algebra_bilinearity_test
  )
  add_test(
    NAME algebra_groups_test
    COMMAND algebra_groups_test
  )
  add_test(
    NAME algebra_fields_test
    COMMAND algebra_fields_test
  )

  add_dependencies(check algebra_bilinearity_test)
  add_dependencies(check algebra_groups_test)
  add_dependencies(check algebra_fields_test)
  add_dependencies(check time_test)

  add_executable(
    multiexp_profile
    EXCLUDE_FROM_ALL

    algebra/scalar_multiplication/multiexp_profile.cpp
  )
  target_link_libraries(
    multiexp_profile

    ${OPENSSL_LIBRARIES}
    ff
  )

  add_dependencies(profile multiexp_profile)
endif()
