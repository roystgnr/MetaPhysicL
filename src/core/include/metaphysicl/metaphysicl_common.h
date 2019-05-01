#ifndef METAPHYSICL_COMMON_H
#define METAPHYSICL_COMMON_H

// The libmesh_dbg_var() macro indicates that an argument to a function
// is used only in debug mode (i.e., when NDEBUG is not defined).
#ifndef NDEBUG
#define metaphysicl_dbg_var(var) var
#else
#define metaphysicl_dbg_var(var)
#endif

#endif // METAPHYSICL_COMMON_H
