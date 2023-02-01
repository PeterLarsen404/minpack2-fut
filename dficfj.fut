-- Generate the x array from nint.
-- This function corresponds to the dficfj function with task = 'XS'
def task_XS (nint : i64) : []f64 =
  let zero   : f64 = 0.0
  let one    : f64 = 1.0
  let two    : f64 = 2.0
  let three  : f64 = 3.0
  let six    : f64 = 6.0
  let twelve : f64 = 12.0

  let h : f64 = 1f64 / (f64.i64 nint)

  -- Compute 8 elements of x, nint times
  let x = flatten (tabulate nint (\i ->
    let xt = f64.i64 i*h
    in [xt*xt*(three-two*xt),
        six*xt*(one-xt),
        six*(one-two*xt),
        -twelve,zero,zero,zero,zero]))
  in x

-- Computes the collocation and continuity equations
--   (4 collocation, 4 continuity).
-- "b" goes from 1 to nint.
def mk_eq [n] (nint_ : i64) (r: f64) (b : i64)
              (rhnfhk : [][][][]f64) (x : [n]f64) : [8]f64 =
  let cpts = 4i64
  let deg = 4i64
  let npi = cpts + deg
  let var = b*npi

  -- Create the values for each k in the collocation equations.
  let col_w_pre_red : [4][5][4]f64 =
    tabulate_3d 4 5 4 (\ k i j ->
      if (j > (4-i)) then
        rhnfhk[k,4+j-i,4+j-i,4-i]*x[var+4+j]
      else
        (rhnfhk[k,j,j,j]*x[var+i+j])+(rhnfhk[k,4+j-i,4+j-i,4-i]*x[var+4+j]))

  -- Reduce the array to get w[deg+1] values.
  let col_w_red : [4][5]f64 =
    tabulate_2d 4 5 (\ i j -> reduce (+) 0 col_w_pre_red[i,j])

  -- Compute the collocation equations.
  let col_eq : [4]f64 =map (\ w -> w[4] - r*(w[1]*w[2]-w[0]*w[3])) col_w_red

  -- Create the values in the collocation equations.
  let con_w_pre_red : [4][4]f64 =
    tabulate_2d 4 4 (\ i j ->
      if (j > (4-i)) then
        rhnfhk[1,0,4+j-i,4-i]*x[var+4+j]
      else
        (rhnfhk[1,0,j,j]*x[var+i+j])+(rhnfhk[1,0,4+j-i,4-i]*x[var+4+j]))

  -- Reduce the array to get w[deg] values.
  let con_w_red : [4]f64 = map (\ e -> reduce (+) 0 e) con_w_pre_red

  -- Compute the continuity equations.
  -- If i is equal to (nint-1) compute the last to elements instead of four.
  let con_eq : [4]f64 =
    if (b < nint_) then
      tabulate 4 (\i -> x[var+cpts+deg+i]-con_w_red[i])
    else
      [con_w_red[0]-1f64,con_w_red[1],0f64,0f64]

  -- Concatenate the equations. (Unsure if better way exists).
  in tabulate 8 (\i -> if i < 4 then col_eq[i] else con_eq[i-4]) :> [8]f64


-- dficfj function.
def dficfj (nint : i64) (x : []f64) (r : f64) : []f64 =
  let cpts : i64 = 4
  let deg  : i64 = 4
  let npi : i64 = cpts + deg
  let dim : i64 = deg + cpts - 1
  let dim_ : i64 = dim+1
  let deg_ : i64 = deg+1

  let h = 1f64 / (f64.i64 nint)

  -- Predefined rho array
  let rho : [4]f64 = [
    0.0694318413734436035,
    0.330009490251541138,
    0.669990539550781250,
    0.930568158626556396
  ]

  -- Store every value of hm in a an array.
  let hms = tabulate deg_ (\i -> h**(f64.i64 i))

  -- Store every value of rhoijh in a 3d array. 5(#hms) x 4(#rho) x 8(#dim_)
  let rhoijhs : [deg_][cpts][dim_]f64 = tabulate_3d deg_ cpts dim_ (\ i j k ->
                                          hms[i]*(rho[j]**(f64.i64 k)))

  -- Store every value of nf in an array.
  let nfs : [dim_]f64 = [1,1,2,2*3,2*3*4,2*3*4*5,2*3*4*5*6,2*3*4*5*6*7]
    :> [dim_]f64

  -- Create the 4d rhnfhk array:
  let rhnfhk : [deg_][cpts][dim_][dim_]f64 = tabulate deg_ (\m ->
    tabulate_3d cpts dim_ dim_ (\ i j k -> rhoijhs[m,i,j]/nfs[k]))

  -- Rearrange rhnfhk to make it fit the one in fortran. (Costly)
  let rhnfhk : [cpts][dim_][dim_][deg_]f64 = tabulate cpts (\i ->
      tabulate_3d dim_ dim_ deg_ (\ j k m -> rhnfhk[m,i,j,k]))

  -- Compute fvec
  let fvec =
    flatten (map (\ i -> mk_eq (nint-1) r i rhnfhk x) (iota nint))

  -- Prepend x[0], x[1] and remove last two elements.
  let fvec = tabulate (npi*nint) (\ i -> if i < 2 then x[i] else fvec[i-2])

  in fvec

-- dficfj test nint 1
-- ==
-- compiled input @ ./dataset/dficfj_1.in
-- output @ ./dataset/dficfj_1.out

-- dficfj test nint 2
-- ==
-- compiled input @ ./dataset/dficfj_2.in
-- output @ ./dataset/dficfj_2.out

-- dficfj test nint 10
-- ==
-- compiled input @ ./dataset/dficfj_10.in
-- output @ ./dataset/dficfj_10.out

-- dficfj test nint 100
-- ==
-- compiled input @ ./dataset/dficfj_100.in
-- output @ ./dataset/dficfj_100.out

-- dficfj test nint 1000
-- ==
-- compiled input @ ./dataset/dficfj_1000.in
-- output @ ./dataset/dficfj_1000.out

-- dficfj test nint 10000
-- ==
-- compiled input @ ./dataset/dficfj_10000.in
-- output @ ./dataset/dficfj_10000.out

-- dficfj test nint 100000.
-- ==
-- compiled input @ ./dataset/dficfj_100000.in
-- output @ ./dataset/dficfj_100000.out


-- dficfj test nint 1000000.
-- ==
-- compiled input @ ./dataset/dficfj_1000000.in
-- output @ ./dataset/dficfj_1000000.out


let main (nint : i64) (x : []f64) (r : f64)  =
  dficfj nint x r
