translate = executable('NEXUStranslate',
                       ['translateNuc.cpp'],
                       dependencies: ncl_dep,
                       install: false)

test('roundTripTest_Translate', python_prog,
  args: [
    test_script,
    '-o',
    translate,
    test_dir / 'TranslateIn',
    test_dir / 'TranslateOut'
  ],
  is_parallel: false
)
