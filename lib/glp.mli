module Types :
  sig
   type problem = unit Ctypes.ptr
   val problem : problem Ctypes.typ
   type direction = Min | Max
   val of_direction : direction -> int
   type bound =
      Upper of float
    | Lower of float
    | Free
    | Double of float * float
    | Fixed of float
   val of_bound : bound -> int * float * float
  type kind = GLP_CV | GLP_IV | GLP_BV
  val of_kind : kind -> int
  type stat = 
    | GLP_BS  
    | GLP_NL of float 
    | GLP_NU of float
    | GLP_NF 
    | GLP_NS of float 
  
  val of_stat : stat -> int * float * float 
  type klass = GLP_RF_GMI | GLP_RF_MIR | GLP_RF_COV | GLP_RF_CL
  val of_klass : klass -> int
  type sel = GLP_DN_BRNCH | GLP_UP_BRNCH | GLP_NO_BRNCH
  val of_sel : sel -> int  
  type fmt = GLP_MPS_DECK | GLP_MPS_FILE
  val of_fmt : fmt -> int
  
  module Smcp :
  sig
    type t
    val t : t Ctypes.structure Ctypes.typ
    val msg_lev : (int, t Ctypes.structure) Ctypes.field
    val meth : (int, t Ctypes.structure) Ctypes.field
    val pricing : (int, t Ctypes.structure) Ctypes.field
    val r_test : (int, t Ctypes.structure) Ctypes.field
    val tol_bnd : (float, t Ctypes.structure) Ctypes.field
    val tol_dj : (float, t Ctypes.structure) Ctypes.field
    val tol_piv : (float, t Ctypes.structure) Ctypes.field
    val obj_ll : (float, t Ctypes.structure) Ctypes.field
    val obj_ul : (float, t Ctypes.structure) Ctypes.field
    val it_lim : (int, t Ctypes.structure) Ctypes.field
    val tm_lim : (int, t Ctypes.structure) Ctypes.field
    val out_frq : (int, t Ctypes.structure) Ctypes.field
    val out_dly : (int, t Ctypes.structure) Ctypes.field
    val presolve : (int, t Ctypes.structure) Ctypes.field
    val foo_bar :
      (float Ctypes_static.carray, t Ctypes.structure) Ctypes.field
  end
  module Iptcp :
  sig
    type i 
    val i : i Ctypes.structure Ctypes.typ
    val msg_lev : (int, i Ctypes.structure) Ctypes.field 
    val ord_alg : (int, i Ctypes.structure) Ctypes.field 
    val foo_bar : 
      (float Ctypes_static.carray, i Ctypes.structure) Ctypes.field
  end

  module Iocp :
  sig
  type c
    val c : c Ctypes.structure Ctypes.typ
    val msg_lev : (int, c Ctypes.structure) Ctypes.field 
    val br_tech : (int, c Ctypes.structure) Ctypes.field 
    val bt_tech : (int, c Ctypes.structure) Ctypes.field 
    val tol_int : (float, c Ctypes.structure) Ctypes.field 
    val tol_obj : (float, c Ctypes.structure) Ctypes.field 
    val tm_lim  : (int, c Ctypes.structure) Ctypes.field 
    val out_frq : (int, c Ctypes.structure) Ctypes.field 
    val out_dly : (int, c Ctypes.structure) Ctypes.field 
    val cb_func : (unit, c Ctypes.structure) Ctypes.field  
    val cb_info : (unit, c Ctypes.structure) Ctypes.field 
    val cb_size : (int, c Ctypes.structure) Ctypes.field
    val pp_tech : (int, c Ctypes.structure) Ctypes.field
    val mip_gap : (float, c Ctypes.structure) Ctypes.field 
    val mir_cuts : (int, c Ctypes.structure) Ctypes.field
    val gmi_cuts :  (int, c Ctypes.structure) Ctypes.field
    val cov_cuts : (int, c Ctypes.structure) Ctypes.field
    val clq_cuts : (int, c Ctypes.structure) Ctypes.field
    val presolve : (int, c Ctypes.structure) Ctypes.field
    val binarize : (int, c Ctypes.structure) Ctypes.field
    val fp_heur : (int, c Ctypes.structure) Ctypes.field
    val ps_heur : (int, c Ctypes.structure) Ctypes.field
    val ps_tm_lim : (int, c Ctypes.structure) Ctypes.field
    val sr_heur : (int, c Ctypes.structure) Ctypes.field
    val use_sol : (int, c Ctypes.structure) Ctypes.field
    val save_sol : (char, c Ctypes.structure) Ctypes.field
    val alien : (int, c Ctypes.structure) Ctypes.field
     val foo_bar :
      (float Ctypes_static.carray, c Ctypes.structure) Ctypes.field
  end
   module Bfcp :
   sig
   type b
    val b : b Ctypes.structure Ctypes.typ
    val msg_lev : (int, b Ctypes.structure) Ctypes.field
    (*val  type : (int, b Ctypes.structure) Ctypes.field*)
    val lu_size : (int, b Ctypes.structure) Ctypes.field
    val piv_tol : (float, b Ctypes.structure) Ctypes.field
    val piv_lim : (int, b Ctypes.structure) Ctypes.field
    val suhl : (int, b Ctypes.structure) Ctypes.field
    val eps_tol : (float, b Ctypes.structure) Ctypes.field
    val ax_gro : (float, b Ctypes.structure) Ctypes.field
    val nfs_max : (int, b Ctypes.structure) Ctypes.field
    val upd_tol : (float, b Ctypes.structure) Ctypes.field
    val nrs_max : (int, b Ctypes.structure) Ctypes.field
    val rs_size : (int, b Ctypes.structure) Ctypes.field
    val foo_bar : (float Ctypes_static.carray, b Ctypes.structure) Ctypes.field
  
  end 
  module Tree :
  sig 
  type tr
  val tr : tr Ctypes.structure Ctypes.typ
  end 
  
   module Attr :
   sig 
   type at
    val at :  at Ctypes.structure Ctypes.typ
    val level : (int, at Ctypes.structure) Ctypes.field
    val origin : (int, at Ctypes.structure) Ctypes.field
    val klass :   (int, at Ctypes.structure) Ctypes.field
    val foo_bar : (float Ctypes_static.carray, at Ctypes.structure) Ctypes.field
  end
  
  module Mpscp :
  sig
  type mp
   val mp : mp Ctypes.structure Ctypes.typ
   val blank :(int, mp Ctypes.structure) Ctypes.field
   val obj_name : (char, mp Ctypes.structure) Ctypes.field
   val  tol_mps : (float, mp Ctypes.structure) Ctypes.field
   val  foo_bar : (float Ctypes_static.carray, mp Ctypes.structure) Ctypes.field
   
   end   


  module Cpxcp : 
  sig 
 type cp
  val cp : cp  Ctypes.structure Ctypes.typ
  val foo_bar : (float Ctypes_static.carray, cp Ctypes.structure) Ctypes.field
  end 
  
 module Tran :
 sig
 type tr
 val  tr : tr Ctypes.structure Ctypes.typ
 end  
   
end







val create_problem : unit -> Types.problem


val set_problem_name : Types.problem -> string -> unit
val set_objective_name : Types.problem -> string -> unit
val set_objective_direction : Types.problem -> Types.direction -> unit
val add_rows :  Types.problem -> int -> int
val add_cols :  Types.problem -> int -> int
val set_row_name : Types.problem -> int -> char -> unit
val set_col_name :  Types.problem -> int -> char -> unit
val set_row_bnds : Types.problem -> int -> Types.bound -> unit
val set_col_bnds : Types.problem -> int -> Types.bound -> unit
val set_obj_coef : Types.problem -> int -> float -> unit
val set_mat_row : Types.problem -> int -> int -> int -> float -> unit
val set_mat_col : Types.problem -> int -> int -> int -> float -> unit
val load_matrix : Types.problem -> int -> int -> int -> float -> unit
val check_dup : int -> int -> int -> int -> int -> int
val del_rows : Types.problem -> int -> int -> unit
val del_cols : Types.problem -> int -> int -> unit
val copy_prob : Types.problem -> Types.problem -> int -> unit
val erase_prob : Types.problem -> unit
val delete_prob : Types.problem -> unit
val  get_prob_name : Types.problem -> char 
val get_obj_name : Types.problem -> char
val get_obj_dir : Types.problem -> int 
val get_num_rows : Types.problem -> int 
val get_num_col : Types.problem -> int
val get_row_name : Types.problem -> int -> char
val get_col_name : Types.problem -> int -> char
val get_row_type : Types.problem -> int -> int  
val glp_get_row_lb : Types.problem -> int -> float 
val get_row_ub : Types.problem -> int -> float 
val get_col_type : Types.problem -> int -> int
val get_col_lb : Types.problem -> int -> float
val get_col_ub : Types.problem -> int -> float
val get_obj_coef : Types.problem -> int -> float
val get_num_nz : Types.problem -> int
val get_mat_row : Types.problem -> int -> int -> float -> int 
val get_mat_col : Types.problem -> int -> int -> float -> int 
val create_index : Types.problem -> unit
val find_row : Types.problem -> char -> int 
val find_col : Types.problem -> char -> int
val delete_index : Types.problem -> unit
val set_rii : Types.problem -> int -> float -> unit
val set_sjj : Types.problem -> int -> float -> unit
val get_rii : Types.problem -> int -> float 
val get_sjj : Types.problem -> int -> float
val scale_prob : Types.problem -> int -> unit
val unscale_prob : Types.problem -> unit
val std_basis : Types.problem -> unit
val adv_basis : Types.problem -> int -> unit
val cpx_basis : Types.problem -> unit
val simplex : Types.problem -> Types.Smcp.t -> int 
val init_smcp : Types.Smcp.t -> unit
val get_status : Types.problem -> int 
val get_prim_stat : Types.problem -> int
val get_dual_stat : Types.problem -> int
val get_obj_val : Types.problem -> float 
val get_row_stat : Types.problem -> int -> int
val get_row_prim : Types.problem -> int -> float
val get_row_dual  : Types.problem -> int -> float
val get_col_stat : Types.problem -> int -> int
val get_col_prim : Types.problem -> int -> float
val get_col_dual :  Types.problem -> int -> float
val get_unbnd_ray :  Types.problem -> int
val interior : Types.problem -> Types.Iptcp.i -> int
val init_iptcp : Types.Iptcp.i -> unit
val ipt_status : Types.problem -> int
val ipt_row_prim : Types.problem -> int -> float
val ipt_row_dual : Types.problem -> int -> float
val ipt_col_prim : Types.problem -> int -> float
val ipt_col_dual : Types.problem -> int -> float
val set_col_kind : Types.problem ->  int -> Types.direction -> unit
val get_col_kind : Types.problem -> int -> int
val get_num_int : Types.problem -> int 
val get_num_bin : Types.problem -> int  
val set_row_stat : Types.problem -> int -> Types.stat -> unit
val set_col_stat : Types.problem -> int -> Types.stat -> unit
val intopt : Types.problem -> Types.Iocp.c -> int
val init_iocp : Types.Iocp.c -> unit
val mip_status : Types.problem -> int
val mip_obj_val : Types.problem -> float
val mip_row_val : Types.problem -> int -> float
val mip_col_val :  Types.problem -> int -> float
val check_kkt : Types.problem -> int -> int -> float -> int -> float -> int -> unit
val print_sol : Types.problem -> char -> int
val read_sol : Types.problem -> char -> int 
val write_sol : Types.problem -> char -> int 
val print_ranges : Types.problem -> int -> int -> int -> char -> int
val print_ipt : Types.problem -> char -> int 
val read_ipt : Types.problem -> char -> int 
val write_ipt : Types.problem -> char -> int
val print_mip : Types.problem -> char -> int
val read_mip :  Types.problem -> char -> int
val write_mip : Types.problem -> char -> int
val bf_exists : Types.problem -> int
val factorize : Types.problem -> int
val bf_updated : Types.problem -> int
val get_bfcp : Types.problem -> Types.Bfcp.b -> unit
val set_bfcp :  Types.problem -> Types.Bfcp.b -> unit 
val get_bhead :  Types.problem ->  int -> int 
val get_row_bind : Types.problem ->  int -> int 
val get_col_bind :  Types.problem ->  int -> int
val ftran : Types.problem ->  float -> unit 
val btran : Types.problem ->  float -> unit  
val warm_up : Types.problem -> int
val eval_tab_row : Types.problem -> int -> int -> float -> int 
val eval_tab_col : Types.problem -> int -> int -> float -> int
val transform_row : Types.problem -> int -> int -> float -> int
val transform_col : Types.problem -> int -> int -> float -> int 
val prim_rtest : Types.problem -> int -> int -> float -> int -> float -> int 
val dual_rtest : Types.problem -> int -> int -> float -> int -> float -> int 
val analyze_bound : Types.problem -> int -> float -> int -> float -> int -> unit
val analyze_coef :  Types.problem -> int -> float -> int -> float -> float -> int -> float -> unit
val ios_reason : Types.Tree.tr -> int
val ios_get_prob : Types.Tree.tr -> Types.problem
val ios_tree_size : Types.Tree.tr  -> int -> int -> int -> unit
val ios_curr_node : Types.Tree.tr -> int
val ios_next_node : Types.Tree.tr -> int -> int
val ios_prev_node : Types.Tree.tr -> int -> int
val ios_up_node : Types.Tree.tr -> int -> int
val ios_node_level : Types.Tree.tr -> int -> int
val ios_node_bound : Types.Tree.tr -> int -> float
val ios_best_node : Types.Tree.tr -> int
val ios_mip_gap : Types.Tree.tr -> float
val ios_node_data : Types.Tree.tr -> int -> unit
val ios_row_attr : Types.Tree.tr -> int -> Types.Attr.at -> unit
val ios_pool_size : Types.Tree.tr -> int
(*val ios_add_row : Types.Tree.tr -> char -> int -> int -> int -> int -> float -> int -> float -> int *)
val ios_del_row : Types.Tree.tr -> int -> unit 
val ios_clear_pool : Types.Tree.tr -> unit
val ios_can_branch : Types.Tree.tr -> int -> int
val ios_branch_upon : Types.Tree.tr -> int -> int -> unit
val ios_select_node : Types.Tree.tr -> int -> unit  
val ios_heur_sol : Types.Tree.tr -> float -> int
val ios_terminate : Types.Tree.tr -> unit   
val init_mpscp : Types.Mpscp.mp -> unit
val read_mps : Types.problem -> int -> Types.Mpscp.mp -> char -> int
val write_mps : Types.problem -> int -> Types.Mpscp.mp -> char -> int
val init_cpxcp : Types.Cpxcp.cp -> unit 
val read_lp : Types.problem -> Types.Cpxcp.cp -> char -> int
val write_lp : Types.problem -> Types.Cpxcp.cp -> char -> int
val read_prob : Types.problem -> int -> char -> int
val write_prob : Types.problem -> int -> char -> int
val mpl_read_model : Types.Tran.tr -> char -> int -> int 
val mpl_read_data : Types.Tran.tr -> char -> int -> int 
val mpl_generate : Types.Tran.tr -> char -> int
val mpl_build_prob : Types.Tran.tr -> Types.problem -> unit 






















































