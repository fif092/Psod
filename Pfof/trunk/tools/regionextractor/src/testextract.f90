Program testextract

  Use modextract
  Use modvariables

  Implicit none

  type = 'cuboid'
  filename = '/data/home/roy/Cosmo/Test/pfof_cube_snap_part_data_testsorted'

  center(1) = 0.23
  center(2) = 0.64
  center(3) = 0.51

  dimensions(1) = 0.15
  dimensions(2) = 0.15
  dimensions(3) = 0.15

  Call extract_cuboid()

  Print *,'Centre:',center
  Print *,'Pos? ',size(pos,1), size(pos,2), maxval(pos(1,:)), maxval(pos(2,:)), maxval(pos(3,:)),  minval(pos(1,:)), minval(pos(2,:)), minval(pos(3,:))


End Program testextract
