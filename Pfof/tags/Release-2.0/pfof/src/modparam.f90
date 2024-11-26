Module modparam

  Use modconst
  Character(len=80) :: root, pathinput
  Character(len=80) :: nameinfo, namepart
  Character(len=3)  :: code_index
  Integer(kind=4)   :: Mmin, Mmax, grpsize
  Real(kind=SP)     :: perco
  Logical           :: outcube, star, metal, dofof, readfromcube, dotimings


End Module modparam
