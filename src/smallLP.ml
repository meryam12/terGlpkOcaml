(* résolution du problème 
     maximize 0.6x1 + 0.5 x2
     subject to 
       x1 + 2x2 <= 1
       3x1 + x2 <= 2 
*) 

open Glp

let _ = true
  
let essaie lp = 
  let lp = Glp.create_problem () in
  Glp.set_problem_name lp "Essai";
  Glp.set_objective_direction lp Types.Max;
  Glp.add_rows lp 2 |> ignore;
  Glp.set_row_bnds lp 1 (Types.Upper 1.0);
  Glp.set_row_bnds lp 2 (Types.Upper 2.0);
  Glp.add_cols lp 2 |> ignore ;
  Glp.set_col_name lp 1 "x1";
  Glp.set_col_bnds lp 1 (Types.Lower 0.0);
  Glp.set_obj_coef lp 1 0.6;
  Glp.set_col_name lp 2 "x2";
  Glp.set_col_bnds lp 2 (Types.Lower 0.0);
  Glp.set_obj_coef lp 2 0.5;
  let coefs =
    [ (1,1,1.);
      (1,2,2.);
      (2,1,3.);
      (2,2,1.)
    ]
  in
  Glp.load_matrix lp coefs;
  let params = Glp.Types.Smcp.(Ctypes.allocate t (Ctypes.make t)) in
  Glp.init_smcp params;  
  Glp.simplex lp params |> ignore;
  let z = Glp.get_obj_val lp in
  let x1 = Glp.get_col_prim lp 1 in
  let x2 = Glp.get_col_prim lp 2 in
  Printf.printf "Z = %f\n" z;
  Printf.printf "x1 = %f\n" x1;
  Printf.printf "x2 = %f\n" x2

