#include "ifxc.h"

const char *
ifxc_strerror(ifxc_status status)
{
  switch(status){
  case IFXC_OK:
    return "success";
  case IFXC_E_INVALID_ARGUMENT:
    return "invalid argument";
  case IFXC_E_ALLOCATION:
    return "allocation failure";
  case IFXC_E_UNKNOWN_FEATURE_SET:
    return "unknown feature set";
  case IFXC_E_INVALID_NSPIN:
    return "invalid spin mode";
  case IFXC_E_MISSING_INPUT:
    return "missing input";
  case IFXC_E_UNSUPPORTED_DERIVATIVE:
    return "unsupported derivative";
  case IFXC_E_UNSUPPORTED_VARIABLE:
    return "unsupported variable";
  case IFXC_E_NOT_INITIALIZED:
    return "not initialized";
  case IFXC_E_INTERNAL:
    return "internal error";
  default:
    return "unknown error";
  }
}
