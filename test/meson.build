
# 1. Test: funkyValidIn/funkyValidOut (with -x flag)
test('buildCheck_1_funky', python_prog,
  args: [
    test_script,
    '-x',
    normalizer,
    test_dir / 'funkyValidIn',
    test_dir / 'funkyValidOut'
  ],
  is_parallel: false
)

# 2. Test: ExternalValidIn/ExternalValidOut (with -e flag)
test('buildCheck_2_ExternalValid', python_prog,
  args: [
    test_script,
    external_flag,
    normalizer, # Replaced path string with executable object
    test_dir / 'ExternalValidIn',
    test_dir / 'ExternalValidOut'
  ],
  is_parallel: false
)

# 3. Test: ExternalInvalidIn (with -i and -e flags)
test('buildCheck_3_ExternalInvalid', python_prog,
  args: [
    test_script,
    '-i',
    '-e',
    normalizer, # Replaced path string with executable object
    test_dir / 'ExternalInvalidIn'
  ],
  is_parallel: false
)

# 4. Test: Simple check (with -a flag)
test('buildCheck_4_SimpleA', python_prog,
  args: [
    test_script,
    '-a',
    normalizer # Replaced path string with executable object
  ],
  is_parallel: false
)

# 5. Test: characters.nex
test('buildCheck_5_characters', python_prog,
  args: [
    test_script,
    normalizer, # Replaced path string with executable object
    data_dir / 'characters.nex',
    test_dir / 'data'
  ],
  is_parallel: false
)

# 6. Test: sample.tre
test('buildCheck_6_sample', python_prog,
  args: [
    test_script,
    normalizer, # Replaced path string with executable object
    data_dir / 'sample.tre',
    test_dir / 'data'
  ],
  is_parallel: false
)

# 7. Test: NTSValidIn/NTSValidOut
test('buildCheck_7_NTSValid', python_prog,
  args: [
    test_script,
    normalizer, # Replaced path string with executable object
    test_dir / 'NTSValidIn',
    test_dir / 'NTSValidOut'
  ],
  is_parallel: false
)

# 8. Test: NCLconverter (with -o, -y, and -e flags)
test('buildCheck_8_NCLconverter', python_prog,
  args: [
    test_script,
    '-o',
    '-y',
    '-e',
    converter, # Replaced path string with executable object
    test_dir / '2NexmlIn',
    test_dir / '2NexmlOut'
  ],
  is_parallel: false
)
