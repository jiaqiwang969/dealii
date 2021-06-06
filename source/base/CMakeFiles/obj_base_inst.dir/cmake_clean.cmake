file(REMOVE_RECURSE
  "CMakeFiles/obj_base_inst"
  "bounding_box.inst"
  "data_out_base.inst"
  "function.inst"
  "function_level_set.inst"
  "function_restriction.inst"
  "function_time.inst"
  "function_tools.inst"
  "geometric_utilities.inst"
  "incremental_function.inst"
  "mpi.inst"
  "mpi_noncontiguous_partitioner.inst"
  "mpi_remote_point_evaluation.inst"
  "partitioner.cuda.inst"
  "partitioner.inst"
  "polynomials_rannacher_turek.inst"
  "symbolic_function.inst"
  "symmetric_tensor.inst"
  "tensor_function.inst"
  "tensor_function_parser.inst"
  "time_stepping.inst"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/obj_base_inst.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
