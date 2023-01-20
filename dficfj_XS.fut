-- Generate the x array from nint. This function corresponds to the dficfj function with task = 'XS'
def task_XS_opt (nint : i64) : []f64 =
  let zero   : f64 = 0.0
  let one    : f64 = 1.0
  let two    : f64 = 2.0
  let three  : f64 = 3.0
  let six    : f64 = 6.0
  let twelve : f64 = 12.0

  let h : f64 = 1f64 / (f64.i64 nint)

  -- Exclusive scan.
  let xts = scan (+) zero (replicate (nint-1i64) h)
  let xts : [nint]f64 = tabulate nint (\ i -> if i == 0 then zero else xts[i-1])
  -- Compute 8 elements of x for each xt in xts
  let x = flatten (map (\xt -> [xt*xt*(three-two*xt),six*xt*(one-xt),six*(one-two*xt),-twelve,zero,zero,zero,zero]) xts)
  in x

-- task_XS test nint 2
-- ==
-- compiled input @ ./dataset/dficfj_XS_2.in
-- output @ ./dataset/dficfj_XS_2.out

-- task_XS test nint 10
-- ==
-- compiled input @ ./dataset/dficfj_XS_10.in
-- output @ ./dataset/dficfj_XS_10.out

-- task_XS test nint 100
-- ==
-- compiled input @ ./dataset/dficfj_XS_100.in
-- output @ ./dataset/dficfj_XS_100.out

-- task_XS test nint 1000
-- ==
-- compiled input @ ./dataset/dficfj_XS_1000.in
-- output @ ./dataset/dficfj_XS_1000.out

-- task_XS test nint 10000
-- ==
-- compiled input @ ./dataset/dficfj_XS_10000.in
-- output @ ./dataset/dficfj_XS_10000.out

let main (nint : i64) =
  task_XS_opt nint