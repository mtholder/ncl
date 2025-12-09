
# The external flag constant
external_flag = '-e'

# The dependency object for the NCL library (assumed to be defined outside this file)
# ncl_dep = dependency(...)

# Define all executables
normalizer = executable('NEXUSnormalizer',
                        ['normalizer.cpp'],
                        dependencies: ncl_dep,
                        install: true) # Changed to 'true' based on bin_PROGRAMS

validator = executable('NEXUSvalidator',
                       ['normalizer.cpp'],
                       cpp_args:['-DJUST_VALIDATE_NEXUS'],
                       dependencies: ncl_dep,
                       install: true) # Changed to 'true' based on bin_PROGRAMS

converter = executable('NCLconverter',
                       ['normalizer.cpp', 'converter.cpp', 'us2ml.cpp'],
                       cpp_args:['-DNCL_CONVERTER_APP'],
                       dependencies: ncl_dep,
                       install: true) # Changed to 'true' based on bin_PROGRAMS

inspector = executable('NEXUSinspector',
                        ['normalizer.cpp'],
                        cpp_args:['-DJUST_REPORT_NEXUS'],
                        dependencies: ncl_dep,
                        install: false)

us2ml = executable('NEX_us2ml',
                    ['normalizer.cpp', 'us2ml.cpp'],
                    cpp_args:['-DTO_NEXML_CONVERTER'],
                    dependencies: ncl_dep,
                    install: false)

union = executable('NEXUSunion',
                    ['normalizer.cpp'],
                    cpp_args: '-DMULTIFILE_NEXUS_READER',
                    dependencies: ncl_dep,
                    install: false)

# 1. Test: characters.nex
test('roundTripTest_1', python_prog,
  args: [
    test_script,
    normalizer, # Meson automatically resolves this to the built executable path
    data_dir / 'characters.nex',
    test_data_dir
  ],
  is_parallel: false
)

# 2. Test: sample.tre
test('roundTripTest_2', python_prog,
  args: [
    test_script,
    normalizer,
    data_dir / 'sample.tre',
    test_data_dir
  ],
  is_parallel: false
)

# 3. Test: NTSValidIn/NTSValidOut
test('roundTripTest_3', python_prog,
  args: [
    test_script,
    normalizer,
    test_dir / 'NTSValidIn',
    test_dir / 'NTSValidOut'
  ],
  is_parallel: false
)

# 4. Test: funkyValidIn/funkyValidOut (with -x flag)
test('roundTripTest_4', python_prog,
  args: [
    test_script,
    '-x',
    normalizer,
    test_dir / 'funkyValidIn',
    test_dir / 'funkyValidOut'
  ],
  is_parallel: false
)

# 5. Test: ExternalValidIn/ExternalValidOut (with EXTERNAL_FLAG)
test('roundTripTest_5', python_prog,
  args: [
    test_script,
    external_flag, # Passes the literal string '-e'
    normalizer,
    test_dir / 'ExternalValidIn',
    test_dir / 'ExternalValidOut'
  ],
  is_parallel: false
)

# 6. Test: ExternalInvalidIn (with -i and -e flags)
test('roundTripTest_6', python_prog,
  args: [
    test_script,
    '-i',
    external_flag,
    normalizer,
    test_dir / 'ExternalInvalidIn'
  ],
  is_parallel: false
)

# 7. Test: Last one (with -a flag)
test('roundTripTest_7', python_prog,
  args: [
    test_script,
    '-a',
    normalizer
  ],
  is_parallel: false
)
