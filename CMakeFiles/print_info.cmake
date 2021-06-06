MESSAGE(
"###
#
#  The following targets are available (invoke by $ make <target>):
#
#    all            - compile the library and all enabled components
#    clean          - remove all generated files
#    install        - install into CMAKE_INSTALL_PREFIX
#
#    info           - print this help message
#    help           - print a list of valid top level targets
#
#    edit_cache     - run ccmake for changing (cached) configuration variables
#                     and reruns the configure and generate phases of CMake
#    rebuild_cache  - rerun the configure and generate phases of CMake
#
#    documentation  - build component 'documentation'
#    examples       - build component 'examples'
#    library        - build component 'library'
#    package        - build binary package
#
#    test           - run a minimal set of tests
#
#    setup_tests    - set up testsuite subprojects
#    prune_tests    - remove all testsuite subprojects
#
#    indent         - indent all headers and source files that changed since the
#                     last commit to master, including untracked ones
#    indent-all     - indent all headers and source files
#
###")