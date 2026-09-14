enable_testing()

add_test(
  NAME program_help
  COMMAND bash ${CMAKE_SOURCE_DIR}/tests/regressions/program_help.sh $<TARGET_FILE:seq2expr>
)

add_test(
  NAME annotator_help
  COMMAND bash ${CMAKE_SOURCE_DIR}/tests/regressions/annotator_help.sh $<TARGET_FILE:seqannot>
)
