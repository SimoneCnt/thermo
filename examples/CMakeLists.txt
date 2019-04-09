######################
# TESTING
######################

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

# Add a binary test.
# Test one of the main targets against one command line.
# Test succedd if binary returned EXIT_SUCCEED.
function(cyg_addtest_bin testname target cmdline)
    execute_process( COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/examples/${testname})
    add_test(
        NAME ${testname}
        COMMAND ${target} ${cmdline} ${ARGN}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/examples/${testname}
    )
    add_dependencies(check ${target})
endfunction()

# Add a compare test.
# After running a cyg_addtest_bin, this test compare the output to a reference
function(cyg_addtest_cmp testname reffile)
    add_test(
        NAME ${testname}_cmp 
        COMMAND ${CMAKE_COMMAND} -E compare_files ${reffile} ${CMAKE_CURRENT_SOURCE_DIR}/${testname}/${reffile}.ref
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/examples/${testname}
    )
endfunction()


cyg_addtest_bin(water thermo -A ${CMAKE_CURRENT_SOURCE_DIR}/water/water.inp -o water.out)
cyg_addtest_cmp(water water.out)

cyg_addtest_bin(diala thermo -A ${CMAKE_CURRENT_SOURCE_DIR}/diala/c7eq.inp -B ${CMAKE_CURRENT_SOURCE_DIR}/diala/c7ax.inp -s 1:1 -o diala.out)
cyg_addtest_cmp(diala diala.out)

cyg_addtest_bin(bhp thermo -A ${CMAKE_CURRENT_SOURCE_DIR}/bhp/bhp1.inp -B ${CMAKE_CURRENT_SOURCE_DIR}/bhp/bhp2.inp -s 1:1 -o bhp.out)
cyg_addtest_cmp(bhp bhp.out)

cyg_addtest_bin(conv thermo -A ${CMAKE_CURRENT_SOURCE_DIR}/conv/rig.inp -B ${CMAKE_CURRENT_SOURCE_DIR}/conv/pps.inp -s 1:1 -o conv.out)
cyg_addtest_cmp(conv conv.out)

cyg_addtest_bin(insulin thermo -A ${CMAKE_CURRENT_SOURCE_DIR}/insulin/monomer.inp -B ${CMAKE_CURRENT_SOURCE_DIR}/insulin/dimer.inp -s 2:1 -o insulin.out)
cyg_addtest_cmp(insulin insulin.out)

cyg_addtest_bin(ethane thermo -A ${CMAKE_CURRENT_SOURCE_DIR}/ethane/ethane.inp -o ethane.out)
cyg_addtest_cmp(ethane ethane.out)

cyg_addtest_bin(benzene_dimer thermo -A ${CMAKE_CURRENT_SOURCE_DIR}/benzene_dimer/benzene_dimer.inp -o benzene_dimer.out)
cyg_addtest_cmp(benzene_dimer benzene_dimer.out)

