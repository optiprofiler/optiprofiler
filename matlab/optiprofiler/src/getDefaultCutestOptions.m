function cutest_options = getDefaultCutestOptions()
    
    cutest_options.(CutestOptionKey.PROBLEM_TYPE.value) = 'u';
    cutest_options.(CutestOptionKey.MINDIM.value) = 1;
    cutest_options.(CutestOptionKey.MAXDIM.value) = 5;
    cutest_options.(CutestOptionKey.MINCON.value) = 0;
    cutest_options.(CutestOptionKey.MAXCON.value) = 10;
    cutest_options.(CutestOptionKey.EXCLUDELIST.value) = {};
end