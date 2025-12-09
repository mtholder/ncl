ncltest = executable('ncltest',
                     ['ncltest.cpp'],
                     dependencies: ncl_dep,
                     install: false)

test('roundTripTest_ncl', python_prog,
     args: [
       test_script,
       ncltest,
       test_dir / 'OldValidIn',
       test_dir / 'OldValidOut'
     ],
     is_parallel: false
)
