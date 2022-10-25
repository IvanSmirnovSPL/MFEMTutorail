#mfem = dependency('mfem', method : 'cmake')
mfem = [cxx.find_library('mfem'), cxx.find_library('umfpack'), cxx.find_library('klu')]
executable('tutorialRun', 'main.cpp', dependencies : mfem)