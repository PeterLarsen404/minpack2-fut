import "dficfj"

-- testfunction for AD
def dficfj_test_ad [n] (x : [n]f64) : [n]f64 =
  let nint = n/8i64
  in dficfj nint x 1f64 :> [n]f64

-- futhark-AD test nint 1
-- ==
-- entry: fwd_J rev_J
-- compiled input @ ./dataset/dficfj_XS_1.out
-- output @ ./dataset/dficfj_jac_1.out

-- futhark-AD test nint 2
-- ==
-- entry: fwd_J rev_J
-- compiled input @ ./dataset/dficfj_XS_2.out
-- auto output

-- futhark-AD test nint 10 (Only for benchmark)
-- ==
-- entry: fwd_J rev_J
-- compiled input @ ./dataset/dficfj_XS_10.out
-- auto output

-- futhark-AD test nint 100 (Only for benchmark)
-- ==
-- entry: fwd_J rev_J
-- compiled input @ ./dataset/dficfj_XS_100.out
-- auto output

-- futhark-AD test nint 1000 (Only for benchmark)
-- ==
-- entry: fwd_J rev_J
-- compiled input @ ./dataset/dficfj_XS_1000.out
-- auto output

entry fwd_J [n] (x : [n]f64)  =
  flatten (tabulate n (\ i ->
    jvp dficfj_test_ad x (replicate n 0 with [i] = 1)))

entry rev_J [n] (x : [n]f64)  =
  flatten (transpose (tabulate n (\ i ->
    vjp dficfj_test_ad x (replicate n 0 with [i] = 1))))
