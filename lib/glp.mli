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
    type klass = GLP_RF_GMI | GLP_RF_MIR | GLP_RF_COV | GLP_RF_CLQ
    val of_klass : klass -> int
    type sel = GLP_DN_BRNCH | GLP_UP_BRNCH | GLP_NO_BRNCH
    val of_sel : sel -> int
    type stat = GLP_BS | GLP_NL | GLP_NU | GLP_NF | GLP_NS
    val of_stat : stat -> int
    type fmt = GLP_MPS_DECK | GLP_MPS_FILE
    val of_fmt : fmt -> int
    type solution = GLP_SOL | GLP_IPT | GLP_MIP
    val of_solution : solution -> int
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
        val set_message_level : 'a -> 'b -> 'c
      end
    module Iptcp :
      sig
        type t
        val t : t Ctypes.structure Ctypes.typ
        val msg_lev : (int, t Ctypes.structure) Ctypes.field
        val ord_alg : (int, t Ctypes.structure) Ctypes.field
        val foo_bar :
          (float Ctypes_static.carray, t Ctypes.structure) Ctypes.field
      end
    module Iocp :
      sig
        type t
        val t : t Ctypes.structure Ctypes.typ
        val msg_lev : (int, t Ctypes.structure) Ctypes.field
        val br_tech : (int, t Ctypes.structure) Ctypes.field
        val bt_tech : (int, t Ctypes.structure) Ctypes.field
        val tol_int : (float, t Ctypes.structure) Ctypes.field
        val tol_obj : (float, t Ctypes.structure) Ctypes.field
        val tm_lim : (int, t Ctypes.structure) Ctypes.field
        val out_frq : (int, t Ctypes.structure) Ctypes.field
        val out_dly : (int, t Ctypes.structure) Ctypes.field
        val cb_func : (unit, t Ctypes.structure) Ctypes.field
        val cb_size : (int, t Ctypes.structure) Ctypes.field
        val pp_tech : (int, t Ctypes.structure) Ctypes.field
        val mip_gap : (float, t Ctypes.structure) Ctypes.field
        val mir_cuts : (int, t Ctypes.structure) Ctypes.field
        val gmi_cuts : (int, t Ctypes.structure) Ctypes.field
        val cov_cuts : (int, t Ctypes.structure) Ctypes.field
        val clq_cuts : (int, t Ctypes.structure) Ctypes.field
        val presolve : (int, t Ctypes.structure) Ctypes.field
        val binarize : (int, t Ctypes.structure) Ctypes.field
        val fp_heur : (int, t Ctypes.structure) Ctypes.field
        val ps_heur : (int, t Ctypes.structure) Ctypes.field
        val ps_tm_lim : (int, t Ctypes.structure) Ctypes.field
        val sr_heur : (int, t Ctypes.structure) Ctypes.field
        val use_sol : (int, t Ctypes.structure) Ctypes.field
        val ptr : 'a -> (unit, t Ctypes.structure) Ctypes.field
        val alien : (int, t Ctypes.structure) Ctypes.field
        val foo_bar :
          (float Ctypes_static.carray, t Ctypes.structure) Ctypes.field
      end
    module Bfcp :
      sig
        type t
        val t : t Ctypes.structure Ctypes.typ
        val msg_lev : (int, t Ctypes.structure) Ctypes.field
        val lu_size : (int, t Ctypes.structure) Ctypes.field
        val piv_tol : (float, t Ctypes.structure) Ctypes.field
        val piv_lim : (int, t Ctypes.structure) Ctypes.field
        val suhl : (int, t Ctypes.structure) Ctypes.field
        val eps_tol : (float, t Ctypes.structure) Ctypes.field
        val max_gro : (float, t Ctypes.structure) Ctypes.field
        val nfs_max : (int, t Ctypes.structure) Ctypes.field
        val upd_tol : (float, t Ctypes.structure) Ctypes.field
        val nrs_max : (int, t Ctypes.structure) Ctypes.field
        val rs_size : (int, t Ctypes.structure) Ctypes.field
        val foo_bar :
          (float Ctypes_static.carray, t Ctypes.structure) Ctypes.field
      end
    module Tree : sig type t val t : t Ctypes.structure Ctypes.typ end
    module Attr :
      sig
        type t
        val t : t Ctypes.structure Ctypes.typ
        val level : (int, t Ctypes.structure) Ctypes.field
        val origin : (int, t Ctypes.structure) Ctypes.field
        val klass : (int, t Ctypes.structure) Ctypes.field
        val foo_bar :
          (float Ctypes_static.carray, t Ctypes.structure) Ctypes.field
      end
    module Mpscp :
      sig
        type t
        val t : t Ctypes.structure Ctypes.typ
        val blank : (int, t Ctypes.structure) Ctypes.field
        val obj_name : (string, t Ctypes.structure) Ctypes.field
        val tol_mps : (float, t Ctypes.structure) Ctypes.field
        val foo_bar :
          (float Ctypes_static.carray, t Ctypes.structure) Ctypes.field
      end
    module Cpxcp :
      sig
        type t
        val t : t Ctypes.structure Ctypes.typ
        val foo_bar :
          (float Ctypes_static.carray, t Ctypes.structure) Ctypes.field
      end
    module Tran : sig type t val t : t Ctypes.structure Ctypes.typ end
  end
val create_problem : unit -> Types.problem
val set_problem_name : Types.problem -> string -> unit
val set_objective_name : Types.problem -> string -> unit
val set_objective_direction : Types.problem -> Types.direction -> unit
val add_rows : Types.problem -> int -> int
val add_cols : Types.problem -> int -> int
val set_row_name : Types.problem -> int -> string -> unit
val set_col_name : Types.problem -> int -> string -> unit
val set_row_bnds : Types.problem -> int -> Types.bound -> unit
val set_col_bnds : Types.problem -> int -> Types.bound -> unit
val set_obj_coef : Types.problem -> int -> float -> unit
val set_mat_row : Types.problem -> int -> int -> int -> float -> unit
val set_mat_col : Types.problem -> int -> int -> int -> float -> unit
val load_matrix : Types.problem -> (int * int * float) list -> unit
val check_dup : int -> int -> int -> int -> int -> int
val del_rows : Types.problem -> int -> int -> unit
val del_cols : Types.problem -> int -> int -> unit
val copy_prob : Types.problem -> Types.problem -> int -> unit
val erase_prob : Types.problem -> unit
val delete_prob : Types.problem -> unit
val get_prob_name : Types.problem -> string
val get_obj_name : Types.problem -> string
val get_obj_dir : Types.problem -> int
val get_num_rows : Types.problem -> int
val get_num_cols : Types.problem -> int
val get_row_name : Types.problem -> int -> string
val get_col_name : Types.problem -> int -> string
val get_row_type : Types.problem -> int -> int
val get_row_lb : Types.problem -> int -> float
val get_row_ub : Types.problem -> int -> float
val get_col_type : Types.problem -> int -> int
val get_col_lb : Types.problem -> int -> float
val get_col_ub : Types.problem -> int -> float
val get_obj_coef : Types.problem -> int -> float
val get_num_nz : Types.problem -> int
val get_mat_row : Types.problem -> int -> int -> float -> int
val get_mat_col : Types.problem -> int -> int -> float -> int
val create_index : Types.problem -> unit
val find_row : Types.problem -> string -> int
val find_col : Types.problem -> string -> int
val delete_index : Types.problem -> unit
val set_rii : Types.problem -> int -> float -> unit
val set_sjj : Types.problem -> int -> float -> unit
val get_rii : Types.problem -> int -> float
val get_sjj : Types.problem -> int -> float
val scale_prob : Types.problem -> int -> unit
val unscale_prob : Types.problem -> unit
val set_row_stat : Types.problem -> int -> Types.stat -> unit
val set_col_stat : Types.problem -> int -> Types.stat -> unit
val std_basis : Types.problem -> unit
val adv_basis : Types.problem -> int -> unit
val cpx_basis : Types.problem -> unit
val simplex :
  Types.problem -> Types.Smcp.t Ctypes.structure Ctypes_static.ptr -> int
val init_smcp : Types.Smcp.t Ctypes.structure Ctypes_static.ptr -> unit
val get_status : Types.problem -> int
val get_prim_stat : Types.problem -> int
val get_dual_stat : Types.problem -> int
val get_obj_val : Types.problem -> float
val get_row_stat : Types.problem -> int -> int
val get_row_prim : Types.problem -> int -> float
val get_row_dual : Types.problem -> int -> float
val get_col_stat : Types.problem -> int -> int
val get_col_prim : Types.problem -> int -> float
val get_col_dual : Types.problem -> int -> float
val get_unbnd_ray : Types.problem -> int
val interior : Types.problem -> Types.Iptcp.t Ctypes.structure -> int
val init_iptcp : Types.Iptcp.t Ctypes.structure -> unit
val ipt_status : Types.problem -> int
val ipt_row_prim : Types.problem -> int -> float
val ipt_row_dual : Types.problem -> int -> float
val ipt_col_prim : Types.problem -> int -> float
val ipt_col_dual : Types.problem -> int -> float
val set_col_kind : Types.problem -> 'a -> Types.kind -> int -> unit
val get_col_kind : Types.problem -> int -> int
val get_num_int : Types.problem -> int
val get_num_bin : Types.problem -> int
val intopt : Types.problem -> Types.Iocp.t Ctypes.structure -> int
val init_iocp : Types.Iocp.t Ctypes.structure -> unit
val mip_status : Types.problem -> int
val mip_obj_val : Types.problem -> float
val mip_row_val : Types.problem -> int -> float
val mip_col_val : Types.problem -> int -> float
val check_kkt :
  Types.problem -> int -> int -> float -> int -> float -> int -> unit
val print_sol : Types.problem -> string -> int
val read_sol : Types.problem -> string -> int
val write_sol : Types.problem -> string -> int
val print_ranges : Types.problem -> int -> int -> int -> string -> int
val print_ipt : Types.problem -> string -> int
val read_ipt : Types.problem -> string -> int
val write_ipt : Types.problem -> string -> int
val print_mip : Types.problem -> string -> int
val read_mip : Types.problem -> string -> int
val write_mip : Types.problem -> string -> int
val bf_exists : Types.problem -> int
val factorize : Types.problem -> int
val bf_updated : Types.problem -> int
val get_bfcp : Types.problem -> Types.Bfcp.t Ctypes.structure -> unit
val set_bfcp : Types.problem -> Types.Bfcp.t Ctypes.structure -> unit
val get_bhead : Types.problem -> int -> int
val get_row_bind : Types.problem -> int -> int
val get_col_bind : Types.problem -> int -> int
val ftran : Types.problem -> float -> unit
val btran : Types.problem -> float -> unit
val warm_up : Types.problem -> int
val eval_tab_row : Types.problem -> int -> int -> float -> int
val eval_tab_col : Types.problem -> int -> int -> float -> int
val transform_row : Types.problem -> int -> int -> float -> int
val transform_col : Types.problem -> int -> int -> float -> int
val prim_rtest : Types.problem -> int -> int -> float -> int -> float -> int
val dual_rtest : Types.problem -> int -> int -> float -> int -> float -> int
val analyze_bound :
  Types.problem -> int -> float -> int -> float -> int -> unit
val analyze_coef :
  Types.problem ->
  int -> float -> int -> float -> float -> int -> float -> unit
val ios_reason : Types.Tree.t Ctypes.structure -> int
val ios_get_prob : Types.Tree.t Ctypes.structure -> Types.problem
val ios_tree_size :
  Types.Tree.t Ctypes.structure -> int -> int -> int -> unit
val ios_curr_node : Types.Tree.t Ctypes.structure -> int
val ios_next_node : Types.Tree.t Ctypes.structure -> int -> int
val ios_prev_node : Types.Tree.t Ctypes.structure -> int -> int
val ios_up_node : Types.Tree.t Ctypes.structure -> int -> int
val ios_node_level : Types.Tree.t Ctypes.structure -> int -> int
val ios_node_bound : Types.Tree.t Ctypes.structure -> int -> float
val ios_best_node : Types.Tree.t Ctypes.structure -> int
val ios_mip_gap : Types.Tree.t Ctypes.structure -> float
val ios_node_data : Types.Tree.t Ctypes.structure -> int -> unit
val ios_row_attr :
  Types.Tree.t Ctypes.structure -> Types.Attr.t Ctypes.structure -> unit
val ios_pool_size : Types.Tree.t Ctypes.structure -> int
val ios_del_row : Types.Tree.t Ctypes.structure -> int -> unit
val ios_clear_pool : Types.Tree.t Ctypes.structure -> unit
val ios_can_branch : Types.Tree.t Ctypes.structure -> int -> int
val ios_branch_upon :
  Types.Tree.t Ctypes.structure -> Types.sel -> int -> unit
val ios_select_node : Types.Tree.t Ctypes.structure -> int -> unit
val ios_heur_sol : Types.Tree.t Ctypes.structure -> float -> int
val ios_terminate : Types.Tree.t Ctypes.structure -> unit
val init_mpscp : Types.Mpscp.t Ctypes.structure -> unit
val read_mps :
  Types.problem ->
  Types.fmt -> 'a -> Types.Mpscp.t Ctypes.structure -> string -> int
val write_mps :
  Types.problem ->
  Types.fmt -> 'a -> Types.Mpscp.t Ctypes.structure -> string -> int
val init_cpxcp : Types.Cpxcp.t Ctypes.structure -> unit
val read_lp :
  Types.problem -> Types.Cpxcp.t Ctypes.structure -> string -> int
val write_lp :
  Types.problem -> Types.Cpxcp.t Ctypes.structure -> string -> int
val read_prob : Types.problem -> int -> string -> int
val write_prob : Types.problem -> int -> string -> int
val mpl_alloc_wksp : unit -> Types.Tran.t Ctypes.structure Ctypes_static.ptr
val mpl_read_model : Types.Tran.t Ctypes.structure -> string -> int -> int
val mpl_read_data : Types.Tran.t Ctypes.structure -> string -> int
val mpl_generate : Types.Tran.t Ctypes.structure -> string -> int
val mpl_build_prob : Types.Tran.t Ctypes.structure -> Types.problem -> unit
val mpl_free_wksp : Types.Tran.t Ctypes.structure -> unit
val main : int -> string -> int
val read_cnfsat : Types.problem -> string -> int
val check_cnfsat : Types.problem -> int
val write_cnfsat : Types.problem -> string -> int
val minisat1 : Types.problem -> int
val init_env : unit -> int
val version : unit -> string
val free_env : unit -> int
val puts : string -> unit
val term_out : int -> int
val open_tee : string -> int
val close_tee : unit -> int
val assert_ : string -> string -> int -> unit
val alloc : int -> int -> unit Ctypes_static.ptr
val realloc : unit Ctypes_static.ptr -> int -> int -> unit
val free : unit Ctypes_static.ptr -> unit
val mem_limit : int -> unit
