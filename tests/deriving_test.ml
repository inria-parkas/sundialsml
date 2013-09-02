open Pprint

type 'a poly1 = 'a
and  ('a,'b) poly2 = 'a * 'b
and  ('a,'b,'c) poly3 = 'a
and  record1 = { mutable mfield11 : unit;
                 field12: unit; }
and 'a poly1_record2 = { opt_field21 : 'a;
                         mutable opt_mfield22 : 'a array;
                         mutable opt_mfield23 : 'a Lazy.t;
                       }
and ('a,'b) poly2_record3 = { polyrec_field31 : ('a,int) poly2_record3 lazy_t;
                              backref_field32 : ('a,'b) poly2;
                              fwdref_field33 : var1_0
                            }
and var1_0 = OptVar11
and var2_00 = Var21 | Var22
and var3_012 = Var31 | OptVar32 of unit | Var33 of unit * unit
and var4_pair = Var41 of (unit * unit)
and 'a poly1_var5_12 = Var51 of 'a | Var52 of 'a * 'a
and ('a,'b) poly2_var6_123 = Var61 of 'a | Var62 of 'a * 'b
                           | OptVar63 of 'a * 'b * 'b
and ('a,'b,'c) poly3_var7_2 = Var71 of ('a * 'b) * 'c
deriving (pretty ~prefix:Deriving_test
                 ~rename:(Lazy.t to lazy_t)
                 ~optional:(Deriving_test,
                            opt_field21,
                            opt_mfield22,
                            opt_mfield23,
                            OptVar11,     (* <- this should have no effect *)
                            OptVar32,
                            OptVar63))

let test name pp x =
  Format.printf "<%s>@\nRWI:    " name;
  with_read_write_invariance (fun () -> pp Format.std_formatter x);
  Format.printf "@\nNo RWI: ";
  pp Format.std_formatter x;
  Format.printf "@\n@\n"
in
test "poly1" (pp_poly1 pp_unit) ();
test "poly1" (pp_poly1 pp_int) (-1);
test "poly2" (pp_poly2 pp_int pp_int) (-1, -1);
test "poly2" (pp_poly2 pp_int pp_int) (1, 1);
test "poly3" (pp_poly3 pp_unit pp_int pp_int) ();
test "record1" pp_record1 { mfield11 = (); field12 = () };
test "poly1_record2" (pp_poly1_record2 pp_unit)
  { opt_field21 = (); opt_mfield22 = [| (); () |];
    opt_mfield23 = lazy () };
test "poly1_record2" (pp_poly1_record2 pp_int)
  { opt_field21 = -1; opt_mfield22 = [| -1; 0 |];
    opt_mfield23 = lazy (5+3) };
test "poly2_record3" (pp_poly2_record3 pp_unit pp_unit)
  { polyrec_field31 = lazy (invalid_arg "oh-oh");
    backref_field32 = (), ();
    fwdref_field33 = OptVar11 };
test "poly2_record3" (pp_poly2_record3 pp_int pp_bool)
  { polyrec_field31 = Lazy.lazy_from_val
        { polyrec_field31 = lazy (invalid_arg "oh-oh");
          backref_field32 = (1, 2);
          fwdref_field33 = OptVar11 };
    backref_field32 = -1, false;
    fwdref_field33 = OptVar11 };
test "var1_0" pp_var1_0 OptVar11;
test "var2_00" (pp_list pp_var2_00) [Var21; Var22];
test "var3_012" (pp_list pp_var3_012) [Var31; OptVar32 (); Var33 ((), ())];
test "var4_pair" pp_var4_pair (Var41 ((), ()));
test "poly1_var5_12" (pp_list (pp_poly1_var5_12 pp_unit))
  [Var51 (); Var52 ((), ())];
test "poly1_var5_12" (pp_list (pp_poly1_var5_12 pp_int))
  [Var51 0; Var51 (-1); Var52 (0, -1)];
test "poly2_var6_123" (pp_list (pp_poly2_var6_123 pp_unit pp_unit))
  [Var61 (); Var62 ((), ()); OptVar63 ((), (), ())];
test "poly2_var6_123" (pp_list (pp_poly2_var6_123 pp_int pp_bool))
  [Var61 1; Var61 (-1); Var62 (0, true); OptVar63 (0, true, false)];
test "poly3_var7_2" (pp_list (pp_poly3_var7_2 pp_unit pp_unit pp_unit))
  [Var71 (((), ()), ())];
test "poly3_var7_2" (pp_list (pp_poly3_var7_2 pp_int pp_bool pp_char))
  [Var71 ((-1, true), 'a')]
;;
