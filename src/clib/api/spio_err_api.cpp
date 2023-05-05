#include "pio_config.h"
#include "pio.h"
#include "pio_internal.h"
#include "pio_api_impl.h"

/* ========== Error handling APIs =========== */
int PIOc_strerror(int pioerr, char *errmsg, size_t errmsg_sz)
{
  return PIOc_strerror_impl(pioerr, errmsg, errmsg_sz);
}

int PIOc_set_log_level(int level)
{
  return PIOc_set_log_level_impl(level);
}


