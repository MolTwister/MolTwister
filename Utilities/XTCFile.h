#if __has_include(<gromacs/fileio/filetypes.h>) && __has_include(<gromacs/fileio/trxio.h>)

#include <gromacs/fileio/filetypes.h>
#include <gromacs/fileio/trxio.h>

// In some cases the standard gromacs include file for XTC files is not included
// in the gromacs install. Hence, we need to check if it is there. If not create our own.
#if __has_include(<gromacs/fileio/xtcio.h>)
#include <gromacs/fileio/xtcio.h>
#else

// All functions return 1 if successful, 0 otherwise
extern struct t_fileio* open_xtc(const char* filename, const char* mode);   // Open a file for xdr I/O
extern void close_xtc(struct t_fileio* fio);                                // Close the file for xdr I/O */
extern int read_first_xtc(struct t_fileio* fio,
                   int*             natoms,
                   int64_t*         step,
                   real*            time,
                   matrix           box,
                   rvec**           x,
                   real*            prec,
                   gmx_bool*        bOK);                                   // Open xtc file, read xtc file first time, allocate memory for x
extern int read_next_xtc(struct t_fileio* fio,
                   int              natoms,
                   int64_t*         step,
                   real*            time,
                   matrix           box,
                   rvec*            x,
                   real*            prec,
                   gmx_bool*        bOK);                                   // Read subsequent frames
extern int write_xtc(struct t_fileio* fio,
                   int              natoms,
                   int64_t          step,
                   real             time,
                   const rvec*      box,
                   const rvec*      x,
                   real             prec);                                  // Write a frame to xtc file

#endif

#define GROMACS_XTC_HEADERS_INSTALLED 1

#endif
