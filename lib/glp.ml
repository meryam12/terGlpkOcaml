open Ctypes
open Foreign

(* For now, let's just put everything in this file. I separate the
   type definitions (module Type) from the glpk functions binded in
   Ocaml.

   We might want latter to re-organize this module (separate it into
   multiple files). But that's ok: refactoring is easy in functional
   programming language.

   How to compile: from project directory (the one containing
   directories src and lib) do: 
     ocamlbuild -use-ocamlfind glp.cmo 
   then glp.cmo is in directory _build/lib

   How to infer the interfacle glp.mli: from project directory do:
     ocamlbuild -use-ocamlfind glp.inferred.mli
   then 
     cp _build/lib/glp.inferred.mli lib/glp.mli
   to replace the current mli with a new version. 


   The header glpk.h is likely in /usr/include/glpk.h
   All its functions should be imported here (maybe not the last part
   on graphs, which can wait). Macros should be transformed into data
   types (see type direction below).

   See Real World Ocaml chapter 19 for more info on Ctypes.

   DISCUSSION:
   I do not keep the prefix glp for glpk functions and
   values. Everything in this file is part of module Glp (the name of
   the file), so will be called prefixed by [Glp.] from outside this
   module (for instance Glp.Min, Glp.create_problem)

   Also, we should correct the terrible naming of glpk, by using full
   names as long as possible. We can admit some abbreviations like lp
   for linear program. Otherwise no abbreviations, no missing vowels,
   and so on.
*)


module Types =
struct 

  (* A glp_prob is a pointer to an abstract type, 
     we represent it as a pointer to void *)
  type problem = unit ptr

  (* the value [prob] can be used to describe the type of C functions
     in calls to [foreign].
  *)
  let problem : problem typ = ptr void

  (* A type representing the direction of optimization (minimize or
     maximize the lp). We do not want to use the ints 1 or 2 as in C,
     as this is not precise (it allows the use of 3 even when 3 does
     not represent a direction).
  *)
  type direction = Min | Max

  let of_direction = function
    | Min -> 1
    | Max -> 2

  type bound = 
    | Upper of float  
    | Lower of float 
    | Free 
    | Double of float * float 
    | Fixed of float 

  let of_bound = function
    | Upper upper -> (1, 0., upper) 
    | Lower lower -> (2, lower, 0.)
    | Free -> (3,0.,0.)
    | Double (lower,upper) -> (4,lower,upper)
    | Fixed value -> (5, value, value)

  type kind = GLP_CV | GLP_IV | GLP_BV 

  let of_kind = function 
    | GLP_CV -> 1
    | GLP_IV -> 2 
    | GLP_BV -> 3


  type klass = GLP_RF_GMI | GLP_RF_MIR | GLP_RF_COV | GLP_RF_CLQ

  let of_klass = function 
    | GLP_RF_GMI -> 1
    | GLP_RF_MIR -> 2
    | GLP_RF_COV -> 3 
    | GLP_RF_CLQ -> 4


  type sel = GLP_DN_BRNCH | GLP_UP_BRNCH | GLP_NO_BRNCH
  let of_sel = function 

    | GLP_DN_BRNCH -> 1
    | GLP_UP_BRNCH -> 2
    | GLP_NO_BRNCH -> 3


 


(*
GLP_BS
basic variable;
GLP_NL
non-basic variable having active lower bound;
GLP_NU
non-basic variable having active upper bound;
GLP_NF
non-basic free variable;
GLP_NS
non-basic fixed variable *)

  type stat = 
    | GLP_BS  
    | GLP_NL
    | GLP_NU
    | GLP_NF 
    | GLP_NS

  let of_stat = function 
    | GLP_BS -> 1
    | GLP_NL -> 2
    | GLP_NU -> 3
    | GLP_NF -> 4
    | GLP_NS -> 5


  (*The parameter
    fmt
    specifies the MPS format version as follows:
    GLP_MPS_DECK
    fixed (ancient) MPS format;
    GLP_MPS_FILE
    free (modern) MPS format. *)
  type fmt = GLP_MPS_DECK | GLP_MPS_FILE
  let of_fmt = function  
    | GLP_MPS_DECK -> 1
    | GLP_MPS_FILE -> 2




(*
GLP_SOL
basic solution;
GLP_IPT
interior-point solution;
GLP_MIP
mixed integer solution *)

 type solution =
  | GLP_SOL
  | GLP_IPT
  | GLP_MIP

 let of_solution = function
  | GLP_SOL -> 1
  | GLP_IPT -> 2
  | GLP_MIP -> 3













  (* the type [Smcp.t] represents the C type [struct glp_smcp].
     the value [Smcp.t] can be used to describe the type of C functions
     in calls to [foreign]
  *)
  module Smcp = struct
    (* smcp : simplex method control parameters
       DISCUSSION: 
       Maybe we should rename this. One possibility would be to make a 
       module Params with submodules like SimplexMethod, InteriorPoint, ...
       The user could write module Smcp = Params.SimplexMethod if he prefers
       the original naming.

    *)
    type t

    (* note how struct types are defined in Ctypes:
       by giving every field in order with its type
    *)
    let t : t structure typ = structure "glp_smcp"
    (* TODO: rename the fields with decent names *)
    let msg_lev = field t "msg_lev" int 
    let meth = field t "meth" int
    let pricing = field t "pricing" int
    let r_test = field t "r_test" int 
    let tol_bnd = field t "tol_bnd" float
    let tol_dj = field t "tol_dj" float
    let tol_piv = field t "tol_piv" float
    let obj_ll = field t "obj_ll" float
    let obj_ul = field t "obj_ul" float
    let it_lim = field t "it_lim" int
    let tm_lim = field t "tm_lim" int 
    let out_frq = field t "out_frq" int
    let out_dly = field t "int out" int
    let presolve = field t "presolve" int 
    let foo_bar = field t "foo_bar" (array 36 double)
    let () = seal t (* close the type: no more fields *)

    let set_message_level params level =
      assert false (* TODO *)
  end 
  module Iptcp = struct 


    type t
    (* note how struct types are defined in Ctypes:
       by giving every field in order with its type
    *)
    let t : t structure typ = structure "glp_iptcp"
    (* TODO: rename the fields with decent names *)
    let msg_lev = field t "msg_lev" int  
    let ord_alg = field t "ord_alg" int
    let foo_bar = field t "foo_bar" (array 48 double) 


    let () = seal t (* close the type: no more fields *)
  end


  module Iocp = struct

    type t
    let t : t structure typ = structure "glp_iocp"
    let msg_lev = field t "msg_lev" int 
    let br_tech = field t "br_tech"int
    let bt_tech = field t "bt_tech" int
    let tol_int = field t "tol_int" float
    let tol_obj = field t "tol_obj" float
    let tm_lim  = field t "out_frq" int
    let out_frq = field t "out_frq" int
    let out_dly = field t "out_dly" int
    let cb_func = field t "cb_func" void
    let ptr cb_info = field t "cb_info" void
    let cb_size = field t "cb_size" int
    let pp_tech = field t "pp_tech" int
    let mip_gap = field t "mip_gap" float
    let mir_cuts = field t "mir_cuts" int
    let gmi_cuts = field t "gmi_cuts" int
    let cov_cuts = field t "cov_cuts" int
    let clq_cuts = field t "clq_cuts" int
    let presolve = field t "presolve" int
    let binarize = field t "binarize" int
    let fp_heur = field t "fp_heur" int
    let ps_heur = field t "ps_heur" int
    let ps_tm_lim = field t "ps_tm_lim" int
    let sr_heur = field t "sr_heur" int
    let use_sol = field t "use_sol" int
    let ptr save_sol = field t "*save_sol" void
    let alien = field t "alien" int
    let foo_bar = field t "foo_bar" (array 24 double)
    let () = seal t (* close the type: no more fields *)
  end


  
  module Bfcp = struct
    type t
    let t : t structure typ = structure "glp_bfcp" 
    let msg_lev = field t "msg_lev" int
    (* let type = field t "type" int*)
    let lu_size = field t "lu_size" int
    let piv_tol = field t "piv_tol" float
    let piv_lim = field t "piv_lim" int
    let suhl = field t "suhl" int 
    let eps_tol = field t "eps_tol" float
    let max_gro = field t "max_gro" float
    let nfs_max = field t "nfs_max" int
    let upd_tol = field t "upd_tol" float
    let nrs_max = field t "nrs_max" int
    let rs_size = field t "rs_size" int
    let foo_bar = field t "foo_bar" (array 38 double) 
    let () = seal t (* close the type: no more fields *)
  end 
  module Tree = struct 
    type t
    let t : t structure typ = structure "glp_tree"
    let () = seal t (* close the type: no more fields *)
  end


  module Attr = struct
    type t
    let t : t structure typ = structure "glp_attr" 
    let level = field t "level" int
    let origin = field t "origin" int
    let klass = field t "klass" int
    let foo_bar = field t "foo_bar" (array 7 double)
    let () = seal t (* close the type: no more fields *)
  end 

  module Mpscp = struct
    type t
    let t : t structure typ = structure "glp_mpscp"
    let blank = field t "blank" int
    let obj_name = field t "obj_name" string
    let tol_mps = field t "tol_mps" float
    let foo_bar = field t "foo_bar" (array 17 double) 
    let () = seal t (* close the type: no more fields *)
  end     

  module Cpxcp = struct
    type t
    let t : t structure typ = structure "glp_cpxcp" 
    let foo_bar = field t "foo_bar" (array 20 double) 
    let () = seal t (* close the type: no more fields *)
  end 

  module Tran = struct 
    type t
    let t : t structure typ = structure "glp_tran" 
    let () = seal t (* close the type: no more fields *)
  end 




end


open Types


(* glp_prob *glp_create_prob(void) *)
let create_problem =
  foreign "glp_create_prob" (void @-> returning problem)


(* void glp_set_prob_name(glp_prob *P, const char *name); *)
let set_problem_name =
  foreign "glp_set_prob_name"
    (problem @-> string @-> returning void)

(* void glp_set_obj_name(glp_prob *P, const char *name); *)
let set_objective_name =
  foreign "glp_set_obj_name"
    (problem @-> string @-> returning void)


(* void glp_set_obj_dir(glp_prob *P, int dir); *)
let set_objective_direction lp direction =
  let go = 
    foreign "glp_set_obj_dir"
      (problem @-> int @-> returning void)
  in
  go lp (of_direction direction)

(* int glp_add_rows(glp_prob *lp, int nrs); *)
let add_rows =
  foreign "glp_add_rows"
    (problem @-> int @->returning int)


(* int glp_add_cols(glp_prob *lp, int ncs); *)
let add_cols =
  foreign "glp_add_cols"
    (problem @-> int @-> returning int)





(* void glp_set_row_name(glp_prob *lp, int i, const char *name); *)

let set_row_name  =
  foreign "glp_set_row_name"
    (problem @-> int @-> string @-> returning void)


(*void glp_set_col_name(glp_prob *lp, int j, const char *name); *)

let set_col_name =
  foreign "glp_set_col_name" 
    ( problem @-> int @-> string @-> returning void)


(*void glp_set_row_bnds(glp_prob *lp, int i, int type, double lb, double ub);*)
let set_row_bnds lp row bound =
  let go =
    foreign "glp_set_row_bnds" 
      (problem @-> int @-> int @-> float @-> float @-> returning void)
  in
  let (kind, lower_b, upper_b) = of_bound bound in 
  go lp row kind lower_b upper_b 


(*void glp_set_col_bnds(glp_prob *lp, int j, int type, double lb, double ub); *)
let set_col_bnds lp column bound =
  let go =
    foreign "glp_set_col_bnds" 
      (problem @-> int @-> int @-> float @-> float @-> returning void)
  in
  let (kind, lower, upper) = of_bound bound in 
  go lp column kind lower upper


(*void glp_set_obj_coef(glp_prob *lp, int j, double coef);*)
let set_obj_coef =
  foreign "glp_set_obj_coef"
    ( problem @-> int @-> float @-> returning void)


(* void glp_set_mat_row(glp_prob *lp, int i, int len, const int ind[], const double val[]);*)
let set_mat_row  =
  foreign "glp_set_mat_row"
    (problem @-> int @-> int @-> int @-> float @-> returning void)



(*void glp_set_mat_col(glp_prob *lp, int j, int len, const int ind[], const double val[]);*)
let set_mat_col  =
  foreign "glp_set_mat_row"
    (problem @-> int @-> int @-> int @-> float @-> returning void)


(*void glp_load_matrix(glp_prob *lp, int ne, const int ia[], const int ja[], const double ar[]);*)
let load_matrix lp non_zeros =
  let go = 
    foreign "glp_load_matrix"
      (problem @-> int @-> ptr int @-> ptr int @-> ptr float @-> returning void)
  in
  let rows =
    non_zeros
    |> List.map (fun (r,c,v) -> r)
    |> (fun l -> 0 :: l)
    |> CArray.of_list int
    |> CArray.start
  in
  let cols =
    non_zeros
    |> List.map (fun (r,c,v) -> c)
    |> (fun l -> 0 :: l)
    |> CArray.of_list int
    |> CArray.start
  in
  let coefs =
    non_zeros
    |> List.map (fun (r,c,v) -> v)
    |> (fun l -> 0. :: l)
    |> CArray.of_list double
    |> CArray.start
  in 
  go lp (List.length non_zeros) rows cols coefs


(* int glp_check_dup(int m, int n, int ne, const int ia[], const int ja[]) ;*)
let check_dup  =
  foreign "glp_check_dup"
    (int @-> int @-> int @-> int @-> int @-> returning int )



(*void glp_del_rows(glp_prob *lp, int nrs, const int num[]);*)
let del_rows =
  foreign "glp_del_rows"
    (problem @-> int @-> int @-> returning void)


(* void glp_del_cols(glp_prob *lp, int ncs, const int num[]); *)
let del_cols =
  foreign "glp_del_cols"
    (problem @-> int @-> int @-> returning void)



(*void glp_copy_prob(glp_prob *dest, glp_prob *prob, int names);*)
let copy_prob =
  foreign "glp_copy_prob"
    (problem @-> problem @-> int @-> returning void )

(*void glp_erase_prob(glp_prob *lp); *)
let erase_prob =
  foreign "glp_erase_prob"
    (problem @-> returning void )


(*void glp_delete_prob(glp_prob *lp);*)
let delete_prob =
  foreign "glp_delete_prob"
    (problem @-> returning void )





(*const char *glp_get_prob_name(glp_prob *P);*)
let get_prob_name =
  foreign "glp_get_prob_name"
    (problem @-> returning string )


(*const char *glp_get_obj_name(glp_prob *P);*)
let get_obj_name =
  foreign "glp_get_obj_name"
    (problem @-> returning string)



(*int glp_get_obj_dir(glp_prob *P);*)
let get_obj_dir =
  foreign "glp_get_obj_dir"
    (problem @-> returning int )

(*int glp_get_num_rows(glp_prob *P);*)
let get_num_rows =
  foreign "glp_get_num_rows"
    (problem @-> returning int)

(*int glp_get_num_cols(glp_prob *P);*)
let get_num_cols =
  foreign "glp_get_num_cols"
    (problem @-> returning int)


(*const char *glp_get_row_name(glp_prob *P, int i);*)
let get_row_name =
  foreign "glp_get_row_name"
    (problem @-> int @-> returning string)

(*const char *glp_get_col_name(glp_prob *P, int j);*) 
let get_col_name =
  foreign "glp_get_col_name"
    (problem @-> int @-> returning string)

(*int glp_get_row_type(glp_prob *P, int i);*)
let get_row_type =
  foreign "glp_get_row_type"
    (problem @-> int @-> returning int)


(*double glp_get_row_lb(glp_prob *P, int i);*)
let get_row_lb =
  foreign "glp_get_row_lb"
    (problem @-> int @-> returning float)



(*double glp_get_row_ub(glp_prob *P, int i);*)
let get_row_ub = 
  foreign "glp_get_row_ub"
    (problem @-> int @-> returning float)

(*int glp_get_col_type(glp_prob *P, int j);*)
let get_col_type =
  foreign "glp_get_col_type"
    (problem @-> int @-> returning int)

(*double glp_get_col_lb(glp_prob *P, int j);*)
let get_col_lb = 
  foreign "glp_get_col_lb"
    (problem @-> int @-> returning float)


(*double glp_get_col_ub(glp_prob *P, int j);*)
let get_col_ub = 
  foreign "glp_get_col_ub"
    (problem @-> int @-> returning float)


(*double glp_get_obj_coef(glp_prob *P, int j);*)
let get_obj_coef =
  foreign "glp_get_obj_coef"
    (problem @-> int @-> returning float)

(*int glp_get_num_nz(glp_prob *P);*)
let get_num_nz =
  foreign "glp_get_num_nz"
    (problem @-> returning int)

(*int glp_get_mat_row(glp_prob *P, int i, int ind[], double val[]);*)
let get_mat_row =
  foreign "glp_get_mat_row"
    (problem @-> int @-> int @-> float @-> returning int )

(*int glp_get_mat_col(glp_prob *P, int j, int ind[], double val[]);*)
let get_mat_col =
  foreign "glp_get_mat_col"
    (problem @-> int @-> int @-> float @-> returning int )


(*void glp_create_index(glp_prob *P);*)
let create_index =
  foreign "glp_create_index"
    (problem @-> returning void)

(*int glp_find_row(glp_prob *P, const char *name);*)
let find_row =
  foreign "glp_find_row"
    (problem @-> string @-> returning int)

(*int glp_find_col(glp_prob *P, const char *name);*)
let find_col =
  foreign "glp_find_col"
    (problem @-> string @-> returning int)


(*void glp_delete_index(glp_prob *P);*)
let delete_index = 
  foreign "glp_delete_index"
    (problem @-> returning void)

(*void glp_set_rii(glp_prob *P, int i, double rii);*)
let set_rii = 
  foreign "glp_set_rii"
    (problem @-> int @-> float @-> returning void )

(*void glp_set_sjj(glp_prob *P, int j, double sjj);*)
let set_sjj =
  foreign "glp_set_sjj"
    (problem @-> int @-> float @-> returning void )

(*double glp_get_rii(glp_prob *P, int i);*)
let get_rii = 
  foreign "glp_get_rii"
    (problem @-> int @-> returning float )

(*double glp_get_sjj(glp_prob *P, int j);*)
let get_sjj =
  foreign "glp_get_sjj"
    (problem @-> int @-> returning float )

(*void glp_scale_prob(glp_prob *P, int flags);*)
let scale_prob =
  foreign "glp_scale_prob"
    (problem @-> int @-> returning void )

(*void glp_unscale_prob(glp_prob *P);*)
let unscale_prob =
  foreign "glp_unscale_prob"
    (problem @-> returning void)


(*void glp_set_row_stat(glp_prob *P, int i, int stat);*)
let set_row_stat lp row stat =
  let go =
    foreign "glp_set_row_stat"
      (problem @-> int @-> int @-> returning void)
  in 
  let kind = of_stat stat in
  go lp row kind 


(*void glp_set_col_stat(glp_prob *P, int i, int stat);*)
let set_col_stat lp column stat =
  let go =
    foreign "glp_set_col_stat"
      (problem @-> int @-> int @-> returning void)
  in 
  let kind = of_stat stat in
  go lp column kind 

(*void glp_std_basis(glp_prob *P);*)
let std_basis =
  foreign "glp_std_basis"
    (problem @-> returning void)


(*void glp_adv_basis(glp_prob *P, int flags);*)
let adv_basis =
  foreign "glp_adv_basis"
    (problem @-> int @-> returning void )


(*void glp_cpx_basis(glp_prob *P);*)
let cpx_basis = 
  foreign "glp_cpx_basis"
    (problem @-> returning void)

(*int glp_simplex(glp_prob *P, const glp_smcp *parm);*)
let simplex =
  foreign "glp_simplex"
    (problem @-> ptr Smcp.t @-> returning int )

(*void glp_init_smcp(glp_smcp *parm)*)
let init_smcp =
  foreign "glp_init_smcp"
    (ptr Smcp.t @-> returning void )


(*int glp_get_status(glp_prob *P);*)
let get_status =
  foreign "glp_get_status"
    (problem @-> returning int)

(*int glp_get_prim_stat(glp_prob *P);*)
let get_prim_stat =
  foreign "glp_get_prim_stat"
    (problem @-> returning int )

(*int glp_get_dual_stat(glp_prob *P);*)
let get_dual_stat =
  foreign "glp_get_dual_stat"
    (problem @-> returning int)


(*double glp_get_obj_val(glp_prob *P);*)
let get_obj_val =
  foreign "glp_get_obj_val"
    (problem @-> returning float)

(*int glp_get_row_stat(glp_prob *P, int i);*)
let get_row_stat =
  foreign "glp_get_row_stat"
    (problem @-> int @-> returning int )

(*double glp_get_row_prim(glp_prob *P, int i);*)
let get_row_prim =
  foreign "glp_get_row_prim"
    (problem @-> int @-> returning  float)




(*double glp_get_row_dual(glp_prob *P, int i);*)
let get_row_dual =
  foreign "glp_get_row_dual"
    (problem @-> int @-> returning  float)

(*int glp_get_col_stat(glp_prob *P, int j);*)
let get_col_stat =
  foreign "glp_get_col_stat"
    (problem @-> int @-> returning int)

(*double glp_get_col_prim(glp_prob *P, int j);*)
let get_col_prim =
  foreign "glp_get_col_prim"
    (problem @-> int @-> returning float)

(*double glp_get_col_dual(glp_prob *P, int j);*)
let get_col_dual =
  foreign "glp_get_col_dual"
    (problem @-> int @-> returning float)



(*int glp_get_unbnd_ray(glp_prob *P);*)
let get_unbnd_ray =
  foreign "glp_get_unbnd_ray"
    (problem @-> returning int)



(*int glp_interior(glp_prob *P, const glp_iptcp *parm);*)
let interior =
  foreign "glp_interior"
    (problem @-> Iptcp.t @-> returning int)

(*void glp_init_iptcp(glp_iptcp *parm);*)
let init_iptcp =
  foreign "glp_init_iptcp"
    (Iptcp.t @-> returning void)

(*int glp_ipt_status(glp_prob *P);*)
let ipt_status =
  foreign "glp_ipt_status"
    (problem @-> returning int)



(*double glp_ipt_row_prim(glp_prob *P, int i);*)
let ipt_row_prim =
  foreign "glp_ipt_row_prim"
    (problem @-> int @-> returning float)

(*double glp_ipt_row_dual(glp_prob *P, int i);*)
let ipt_row_dual =
  foreign "glp_ipt_row_dual"
    (problem @-> int @-> returning float)

(*double glp_ipt_col_prim(glp_prob *P, int j);*)
let ipt_col_prim =
  foreign "glp_ipt_col_prim"
    (problem @-> int @-> returning float)


(*double glp_ipt_col_dual(glp_prob *P, int j);*)
let ipt_col_dual =
  foreign "glp_ipt_col_dual"
    (problem @-> int @-> returning float)

(*void glp_set_col_kind(glp_prob *P, int j, int kind);*)
let set_col_kind lp j kind =
  let go =
    foreign "glp_set_col_kind"
      (problem@-> int @-> int @-> returning void )
  in
  go lp (of_kind kind)



(*int glp_get_col_kind(glp_prob *P, int j);*)
let get_col_kind =
  foreign "glp_get_col_kind"
    (problem @-> int @-> returning int)

(*int glp_get_num_int(glp_prob *P);*)
let get_num_int =
  foreign "glp_get_num_int"
    (problem @-> returning int)

(*int glp_get_num_bin(glp_prob *P);*)
let get_num_bin =
  foreign "glp_get_num_bin"
    (problem @-> returning int) 








(*int glp_intopt(glp_prob *P, const glp_iocp *parm);*)
let intopt =
  foreign "glp_intopt" 
    (problem @-> Iocp.t @-> returning int )

(*void glp_init_iocp(glp_iocp *parm);*)
let init_iocp =
  foreign "glp_init_iocp"
    (Iocp.t @-> returning void)

(*int glp_mip_status(glp_prob *P);*)
let mip_status =
  foreign "glp_mip_status"
    (problem @-> returning int)

(*double glp_mip_obj_val(glp_prob *P);*)
let mip_obj_val =
  foreign "glp_mip_obj_val"
    (problem @-> returning float)

(*double glp_mip_row_val(glp_prob *P, int i);*)
let mip_row_val =
  foreign "glp_mip_row_val"
    (problem @-> int @-> returning float)

(*double glp_mip_col_val(glp_prob *P, int j);*)
let mip_col_val =
  foreign "glp_mip_col_val"
    (problem @-> int @-> returning float)


(*void glp_check_kkt(glp_prob *P, int sol, int cond, double *ae_max,
      int *ae_ind, double *re_max, int *re_ind);*)
let check_kkt =
  foreign "glp_check_kkt"
    (problem @-> int @-> int @-> float @-> int @-> float @-> int @-> returning void)


(*int glp_print_sol(glp_prob *P, const char *fname);*)
let print_sol =
  foreign "glp_print_sol"
    (problem @-> string @-> returning int)

(*int glp_read_sol(glp_prob *P, const char *fname);*)
let read_sol =
  foreign "glp_read_sol" 
    (problem @-> string @-> returning int)

(*int glp_write_sol(glp_prob *P, const char *fname);*)
let write_sol =
  foreign "glp_write_sol"
    (problem @-> string @-> returning int)

(*int glp_print_ranges(glp_prob *P, int len, const int list[],
      int flags, const char *fname);*)
let print_ranges =
  foreign "glp_print_ranges"
    (problem @-> int @-> int @-> int @-> string @-> returning int )

(*int glp_print_ipt(glp_prob *P, const char *fname);*)
let print_ipt =
  foreign "glp_print_ipt"
    (problem @-> string @-> returning int  )

(*int glp_read_ipt(glp_prob *P, const char *fname);*)
let read_ipt = 
  foreign "glp_read_ipt"
    (problem @-> string @-> returning int)

(*int glp_write_ipt(glp_prob *P, const char *fname);*)
let write_ipt =
  foreign "glp_write_ipt"
    (problem @-> string @-> returning int)


(*int glp_print_mip(glp_prob *P, const char *fname);*)
let print_mip =
  foreign "glp_print_mip"
    (problem @-> string @-> returning int)

(*int glp_read_mip(glp_prob *P, const char *fname);*)
let read_mip =
  foreign "glp_read_mip"
    (problem @-> string @-> returning int)

(*int glp_write_mip(glp_prob *P, const char *fname);*)
let write_mip =
  foreign "glp_write_mip"
    (problem @-> string @-> returning int)




(*int glp_bf_exists(glp_prob *P);*)
let bf_exists =
  foreign "glp_bf_exists"
    (problem @-> returning int)

(*int glp_factorize(glp_prob *P);*)
let factorize =
  foreign "glp_factorize"
    (problem @-> returning int)

(*int glp_bf_updated(glp_prob *P);*)
let bf_updated =
  foreign "glp_bf_updated"
    (problem @-> returning int)

(*void glp_get_bfcp(glp_prob *P, glp_bfcp *parm);*)
let get_bfcp =
  foreign "glp_get_bfcp" 
    (problem @-> Bfcp.t @-> returning void)

(*void glp_set_bfcp(glp_prob *P, const glp_bfcp *parm);*)
let set_bfcp =
  foreign "glp_set_bfcp"
    (problem @-> Bfcp.t @-> returning void )



(*int glp_get_bhead(glp_prob *P, int k);*)
let get_bhead =
  foreign "glp_get_bhead"
    (problem @-> int @-> returning int)


(*int glp_get_row_bind(glp_prob *P, int i);*)
let get_row_bind =
  foreign "glp_get_row_bind"
    (problem @-> int @-> returning int)

(*int glp_get_col_bind(glp_prob *P, int j);*)
let get_col_bind =
  foreign "glp_get_col_bind"
    (problem @-> int @-> returning int)

(*void glp_ftran(glp_prob *P, double x[]);*)
let ftran =
  foreign "glp_ftran"
    (problem @-> float @-> returning void)


(*void glp_btran(glp_prob *P, double x[]);*)
let btran =
  foreign "glp_btran"
    (problem @-> float @-> returning void)

(*int glp_warm_up(glp_prob *P);*)
let warm_up =
  foreign "glp_warm_up"
    (problem @-> returning int)

(*int glp_eval_tab_row(glp_prob *P, int k, int ind[], double val[]);*)
let eval_tab_row =
  foreign "glp_eval_tab_row"
    (problem @-> int @-> int @-> float @-> returning int)

(*int glp_eval_tab_col(glp_prob *P, int k, int ind[], double val[]);*)
let eval_tab_col =
  foreign "glp_eval_tab_col"
    (problem @-> int @-> int @-> float @-> returning int)


(*int glp_transform_row(glp_prob *P, int len, int ind[], double val[]);*)
let transform_row =
  foreign "glp_transform_row"
    (problem @-> int @-> int @-> float @-> returning int)





(*int glp_transform_col(glp_prob *P, int len, int ind[], double val[]);*)
let transform_col =
  foreign "glp_transform_col"
    (problem @-> int @-> int @-> float @-> returning int)

(*int glp_prim_rtest(glp_prob *P, int len, const int ind[],
      const double val[], int dir, double eps);*)
let prim_rtest =
  foreign "glp_prim_rtest"
    (problem @-> int @-> int @-> float @-> int @-> float @-> returning int)

(*int glp_dual_rtest(glp_prob *P, int len, const int ind[],
      const double val[], int dir, double eps);*)
let dual_rtest =
  foreign "glp_dual_rtest"
    (problem @-> int @-> int @-> float @-> int @-> float @-> returning int)

(*void glp_analyze_bound(glp_prob *P, int k, double *value1, int *var1,
      double *value2, int *var2);*)
let analyze_bound =
  foreign "glp_analyze_bound"
    (problem @-> int @-> float @-> int @-> float @-> int @-> returning void)

(*void glp_analyze_coef(glp_prob *P, int k, double *coef1, int *var1,
      double *value1, double *coef2, int *var2, double *value2);*)
let analyze_coef =
  foreign "glp_analyze_coef"
    (problem @-> int @-> float @-> int @-> float @-> float @-> int @-> float @-> returning void)



(*int glp_ios_reason(glp_tree *T);*)
let ios_reason =
  foreign "glp_ios_reason"
    (Tree.t @-> returning int)

(*glp_prob *glp_ios_get_prob(glp_tree *T);*)
let ios_get_prob =
  foreign "glp_ios_get_prob"
    (Tree.t @-> returning problem)




(*void glp_ios_tree_size(glp_tree *T, int *a_cnt, int *n_cnt,
      int *t_cnt);*)
let ios_tree_size = 
  foreign "glp_ios_tree_size"
    (Tree.t @-> int @-> int @-> int @-> returning void)

(*int glp_ios_curr_node(glp_tree *T);*)
let ios_curr_node = 
  foreign "glp_ios_curr_node"
    (Tree.t @-> returning int)

(*int glp_ios_next_node(glp_tree *T, int p);*)
let ios_next_node =
  foreign "glp_ios_next_node"
    (Tree.t @-> int @-> returning int)

(*int glp_ios_prev_node(glp_tree *T, int p);*)
let ios_prev_node = 
  foreign "glp_ios_prev_node"
    (Tree.t @-> int @-> returning int )

(*int glp_ios_up_node(glp_tree *T, int p);*)
let ios_up_node =
  foreign "glp_ios_up_node"
    (Tree.t @-> int @-> returning int) 

(*int glp_ios_node_level(glp_tree *T, int p);*)
let ios_node_level =
  foreign "glp_ios_node_level"
    (Tree.t @-> int @-> returning int)

(*double glp_ios_node_bound(glp_tree *T, int p);*)
let ios_node_bound =
  foreign "glp_ios_node_bound"
    (Tree.t @-> int @-> returning float)

(*int glp_ios_best_node(glp_tree *T);*)
let ios_best_node =
  foreign "glp_ios_best_node"
    (Tree.t @-> returning int)

(*double glp_ios_mip_gap(glp_tree *T);*)
let ios_mip_gap =
  foreign "glp_ios_mip_gap"
    (Tree.t @-> returning float)

(*void *glp_ios_node_data(glp_tree *T, int p);*)
let ios_node_data =
  foreign "glp_ios_node_data"
    (Tree.t @-> int @-> returning void)


(*void glp_ios_row_attr(glp_tree *T, int i, glp_attr *attr);*)
let ios_row_attr = 
  foreign "glp_ios_row_attr"
    (Tree.t @-> Attr.t @-> returning void)


(*int glp_ios_pool_size(glp_tree *T);*)
let ios_pool_size =
  foreign "glp_ios_pool_size"
    (Tree.t @-> returning int)



(* ERREUR
   (*int glp_ios_add_row(glp_tree *T,
      const char *name, int klass, int flags, int len, const int ind[],
      const double val[], int type, double rhs);*)
   let ios_add_row tree klass =
   let go =
   foreign "glp_ios_add_row"
    (Tree.t @-> string @-> int @-> int @-> int @-> int @-> float @-> int @-> float @-> returning int )
   in
   go tree (of_klass klass)
*)

(*void glp_ios_del_row(glp_tree *T, int i);*)
let ios_del_row = 
  foreign "glo_ios_del_row"
    (Tree.t @-> int @-> returning void)

(*void glp_ios_clear_pool(glp_tree *T);*)
let ios_clear_pool =
  foreign "glp_ios_clear_pool"
    (Tree.t @-> returning void)

(*int glp_ios_can_branch(glp_tree *T, int j);*)
let ios_can_branch =
  foreign "glp_ios_can_branch"
    (Tree.t @-> int @-> returning int)

(*void glp_ios_branch_upon(glp_tree *T, int j, int sel);*)
let ios_branch_upon tree sel =
  let go = 
    foreign "glp_ios_branch_upon"
      (Tree.t @-> int @-> int @-> returning void)
  in 
  go tree (of_sel sel)

(*void glp_ios_select_node(glp_tree *T, int p);*)
let ios_select_node = 
  foreign "glp_ios_select_node"
    (Tree.t @-> int @-> returning void)


(*int glp_ios_heur_sol(glp_tree *T, const double x[]);*)
let ios_heur_sol =
  foreign "glp_ios_heur_sol"
    (Tree.t @-> float @-> returning int)

(*void glp_ios_terminate(glp_tree *T);*)
let ios_terminate = 
  foreign "glp_ios_terminate"
    (Tree.t @-> returning void)


(*void glp_init_mpscp(glp_mpscp *parm);*)
let init_mpscp =
  foreign "glp_init_mpscp"
    (Mpscp.t @-> returning void)

(*int glp_read_mps(glp_prob *P, int fmt, const glp_mpscp *parm,
      const char *fname);*)
let read_mps lp fmt mp = 
  let go =
    foreign "glp_read_mps"
      (problem @-> int @-> Mpscp.t @-> string @-> returning int )
  in 
  go lp (of_fmt fmt)




(*
int glp_write_mps(glp_prob *P, int fmt, const glp_mpscp *parm,
      const char *fname);*)
let write_mps lp fmt mp = 
  let go =
    foreign "glp_write_mps"
      (problem @-> int @-> Mpscp.t @-> string @-> returning int )
  in 
  go lp (of_fmt fmt)


(* void glp_init_cpxcp(glp_cpxcp *parm); *)
let init_cpxcp =
  foreign "glp_init_cpxcp"
    (Cpxcp.t @-> returning void)


(*int glp_read_lp(glp_prob *P, const glp_cpxcp *parm, const char *fname);*)
let read_lp =
  foreign "glp_read_lp"
    (problem @-> Cpxcp.t @-> string @-> returning int )

(*int glp_write_lp(glp_prob *P, const glp_cpxcp *parm, const char *fname);*)
let write_lp =
  foreign "glp_write_lp"
    (problem @-> Cpxcp.t @-> string @-> returning int )

(*int glp_read_prob(glp_prob *P, int flags, const char *fname);*)
let read_prob =
  foreign "glp_read_prob"
    (problem @-> int @-> string @-> returning int)


(*int glp_write_prob(glp_prob *P, int flags, const char *fname);*)
let write_prob =
  foreign "glp_write_prob"
    (problem @-> int @-> string @-> returning int)

(* I don't know
   glp_tran *glp_mpl_alloc_wksp(void);
   /* allocate the MathProg translator workspace */*)
let mpl_alloc_wksp =
  foreign "glp_mpl_alloc_wksp"
    (void @-> returning (ptr Tran.t))

(*int glp_mpl_read_model(glp_tran *tran, const char *fname, int skip);*)
let mpl_read_model =
  foreign "glp_mpl_read_model"
    (Tran.t @-> string @-> int @-> returning int)

(*int glp_mpl_read_data(glp_tran *tran, const char *fname);*)
let mpl_read_data =
  foreign "glp_mpl_read_data"
    (Tran.t @-> string @-> returning int)


(*int glp_mpl_generate(glp_tran *tran, const char *fname);*)
let mpl_generate =
  foreign "glpp_mpl_generate"
    (Tran.t @-> string @-> returning int)


(*void glp_mpl_build_prob(glp_tran *tran, glp_prob *prob);*)
let mpl_build_prob =
  foreign "glp_mpl_build_prob"
    (Tran.t @-> problem @-> returning void)


(*int glp_mpl_postsolve(glp_tran *tran, glp_prob *prob, int sol);
let mpl_postsolve lp =
 let go =
   foreign "glp_mpl_postsolve"
     (ptr Tran.t @-> problem @-> int @-> returning int   )
 in 
 go lp (of_solution solution)

*)


(*void glp_mpl_free_wksp(glp_tran *tran);*)
let mpl_free_wksp =
 foreign "glp_mpl_free_wksp"
   (Tran.t  @-> returning void)

(*int glp_main(int argc, const char *argv[]);*)
let main =
 foreign "glp_main"
  (int @-> string @-> returning int)

(*int glp_read_cnfsat(glp_prob *P, const char *fname);*)
let read_cnfsat =
 foreign "glp_read_cnfsat"
  (problem @-> string @-> returning int)

(*int glp_check_cnfsat(glp_prob *P);*)
let check_cnfsat =
 foreign "glp_check_cnfsat"
  (problem @-> returning int )

(*int glp_write_cnfsat(glp_prob *P, const char *fname);*)
let write_cnfsat = 
 foreign "glp_write_cnfsat"
  (problem @-> string @-> returning int)

(*int glp_minisat1(glp_prob *P);*)
let minisat1 =
  foreign "glp_minisat1"
   (problem @-> returning int)



(*int glp_init_env(void);*)
let init_env =
 foreign "glp_init_env"
   (void @-> returning int)


(*const char *glp_version(void);*)
let version =
 foreign "glp_version"
  (void @-> returning string)

(*int glp_free_env(void);*)
let free_env =
 foreign "glp_free_env"
  (void @-> returning int)

(*void glp_puts(const char *s);*)
let puts =
  foreign "glp_puts"
    (string @-> returning void)





(*int glp_term_out(int flag);*)
let term_out =
 foreign "glp_term_out"
   (int @-> returning int)







(*int glp_open_tee(const char *name);*)
let open_tee =
 foreign "glp_open_tee"
  (string @-> returning int)


(*int glp_close_tee(void);*)
let close_tee =
 foreign "glp_close_tee"
  (void @-> returning int)



(*void glp_assert_(const char *expr, const char *file, int line);*)
let assert_ =
 foreign "glp_assert"
  (string @-> string @-> int @-> returning void)





(*void *glp_alloc(int n, int size);*)
let alloc = 
 foreign "glp_alloc"
  (int @-> int @-> returning (ptr void))

(*void *glp_realloc(void *ptr, int n, int size);*)
let realloc =
 foreign "glp_realloc"
  (ptr void @-> int @-> int @-> returning void)

(*void glp_free(void *ptr);*)
let free =
 foreign "glp_free"
   (ptr void @-> returning void)

(*void glp_mem_limit(int limit);*)
let mem_limit =
 foreign "glp_mem_limit"
  (int @-> returning void)





