module Types :
  sig
    type problem = unit Ctypes.ptr
    val problem : problem Ctypes.typ
    type direction = Min | Max
    val of_direction : direction -> int
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
  end
val create_problem : unit -> Types.problem
val set_problem_name : Types.problem -> string -> unit
val set_objective_name : Types.problem -> string -> unit
val set_objective_direction : Types.problem -> Types.direction -> unit
