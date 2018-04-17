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
   then glp.cmo is in directory _build

   How to infer the interfacle glp.mli: from project directory do:
     ocamlbuild -use-ocamlfind glp.inferred.mli
   then 
     cp _build/glp.inferred.mli lib/glp.mli
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
    let tol_bnd = field t "tol_bnd" double
    let tol_dj = field t "tol_dj" double
    let tol_piv = field t "tol_piv" double
    let obj_ll = field t "obj_ll" double
    let obj_ul = field t "obj_ul" double
    let it_lim = field t "it_lim" int
    let tm_lim = field t "tm_lim" int 
    let out_frq = field t "out_frq" int
    let out_dly = field t "int out" int
    let presolve = field t "presolve" int 
    let foo_bar = field t "foo_bar" (array 36 double)
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


