
GDRIVE := /run/user/1000/gvfs/google-drive:host=gmail.com,user=USERNAME/GVfsSharedWithMe/PATHHASH/OPTIONALMOREHASH

OUTDIR := output

${OUTDIR}:
	ln -s ${GDRIVE} $@
