Program testextract

  Use modextract
  Use modvariables

  Implicit none

  type = 'cuboid'
  filename = '/home/roy/Travail/Devel/Cosmologie/test/trunk/pfof_cube_snap_part_data_testtrunk'

  center(1) = 0.23
  center(2) = 0.64
  center(3) = 0.51

  dimensions(1) = 0.10
  dimensions(2) = 0.10
  dimensions(3) = 0.10

  Call extract_cuboid()

  Print *,maxval(pos(1,:)), maxval(pos(2,:)), maxval(pos(3,:))

  Deallocate(pos)

End Program testextract
