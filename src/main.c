#define _WIN32_IE 0x501
#define _WIN32_WINNT 0x501

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <ocidl.h>
#include <olectl.h>
#include <wingdi.h>
#include <windows.h>
#include <shlwapi.h>
#include <commdlg.h>
#include <commctrl.h>
#include "res\resource.h"

#define BTN_SIZE 6
#define LIN_SIZE (BTN_SIZE / 2)
#define CLR_BACK 0x000000
#define DLG_WNDC "#32770"

#define tr(f) ((LONG)(f))
#define ZZ sizeof(BGRA)



typedef struct _CONF* PCONF;
typedef void FLTRFUNC(PCONF, DWORD);

typedef struct _PICT {
    HDC devc;
    HBITMAP hdib;
    HGLOBAL bptr;
    POINT size;
    struct _PICT *prev, *next;
} PICT;

typedef struct _CONF {
    HWND hwid;
    RECT area;
    HBRUSH mark;
    DWORD data;
    FLOAT dprg, cprg, mprg;
    FLTRFUNC *func;
    PCONF next;
} CONF;

typedef struct _FLTR {
    BYTE tdlg, bdef;
    SHORT liml, limh;
    FLTRFUNC *func;
    LPSTR name;
} FLTR;

typedef struct _TCNF {
    LONG dmin, dmax;
    PCONF conf;
} TCNF;

typedef DWORD CALLBACK THRDFUNC(TCNF*);

typedef struct _HIST {
    DWORD npix, rpix[256], gpix[256], bpix[256];
} HIST;

#pragma pack(push, 1)
typedef union _BGRA {
    struct {
        BYTE B, G, R, A;
    };
    DWORD BGRA;
} BGRA;

typedef struct _BMPH {
    WORD bmp;
    DWORD bsz, rsv, off, img, bmw, bmh;
    WORD pln, bpp;
    BYTE tag[24];
} BMPH;
#pragma pack(pop)



HINSTANCE base;
HBRUSH bkBrush;
BOOL flgPaint;
HCURSOR curStd, curMove, curSize;
HWND mainSurface, fwdButton, bwdButton, progBar, idRect = NULL;
POINT dlgMinSz, pbOffset, ibOffset, cpos, ppos;
DWORD numFreeCores;
PICT *root, *curr;
CONF *conf = NULL;



PICT *MakePict();
void AssignDIB(PICT *pict, HBITMAP hbm, LONG w, LONG h);
void DuplicatePict(PICT *pict);
void FreePictTree(PICT **root);



inline void Progress(CONF *vcnf) {
    if (vcnf->hwid) {
        vcnf->dprg += 1.0;
        if (tr(vcnf->dprg * vcnf->mprg) > vcnf->cprg)
            SendMessage(progBar, PBM_SETPOS, vcnf->cprg = vcnf->dprg * vcnf->mprg, 0);
    }
}



void CutPict(CONF *vcnf) {
    if (vcnf->area.left != 0 || vcnf->area.top != 0 || vcnf->area.right != curr->size.x || vcnf->area.bottom != curr->size.y) {
        BYTE *src, *dst, *dpos;
        LONG x, y;

        curr = curr->prev;
        FreePictTree(&curr->next);
        curr->next = MakePict();
        curr->next->prev = curr;
        curr = curr->next;
        curr->size.x = vcnf->area.right - vcnf->area.left;
        curr->size.y = vcnf->area.bottom - vcnf->area.top;
        AssignDIB(curr, curr->prev->hdib, curr->size.x, -curr->size.y);
        dst = curr->bptr + (curr->size.y * curr->size.x) * ZZ;
        src = curr->prev->bptr;
        vcnf->mprg = 100.0/(vcnf->area.bottom - vcnf->area.top);

        for (y = vcnf->area.bottom - 1; y >= vcnf->area.top; y--) {
            dpos = src + y * curr->prev->size.x * ZZ;
            for (x = vcnf->area.right - 1; x >= vcnf->area.left; x--) {
                dst -= ZZ;
                ((BGRA*)dst)->BGRA = ((BGRA*)(dpos + x * ZZ))->BGRA;
            }
            Progress(vcnf);
        }
    }
}



int MessageBoxIcon(HWND hwnd, LPCSTR text, LPCSTR capt, UINT style, LPCSTR icon) {
    MSGBOXPARAMS mbp;

    mbp.cbSize = sizeof(MSGBOXPARAMS);
    mbp.hwndOwner = hwnd;
    mbp.hInstance = base;
    mbp.lpszText = text;
    mbp.lpszCaption = capt;
    mbp.dwStyle = style;
    mbp.lpszIcon = icon;
    mbp.dwContextHelpId = 0;
    mbp.lpfnMsgBoxCallback = NULL;
    mbp.dwLanguageId = LANG_NEUTRAL;

    return MessageBoxIndirect(&mbp);
}



DWORD Parallelize(CONF *vcnf, LONG dmin, LONG dmax, THRDFUNC func) {
    LONG i, d, n = min(numFreeCores, dmax - dmin);
    HANDLE *htrd = malloc(n * sizeof(HANDLE));
    TCNF *tcnf = malloc(n * sizeof(TCNF));
    DWORD retn;

    d = (dmax - dmin) / n;
    for (i = 0; i < n; i++) {
        tcnf[i].conf = vcnf;
        tcnf[i].dmin = dmin + i * d;
        tcnf[i].dmax = (i == n - 1)? dmax : dmin + (i + 1) * d;
        htrd[i] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)func, &tcnf[i], 0, NULL);
    }
    retn = WaitForMultipleObjects(n, htrd, TRUE, INFINITE);
    free(tcnf);
    free(htrd);
    return retn;
}



void GrayWorld(CONF *vcnf, DWORD iact) {
    if (iact) {
        FLOAT CR, CG, CB, X;
        LONGLONG RR = 0, GG = 0, BB = 0;
        LONG x, y, M = 0;
        BYTE *dpos, *src = curr->bptr;

        vcnf->mprg = 50.0/(vcnf->area.bottom - vcnf->area.top);
        for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
            dpos = src + (vcnf->area.left + (y * curr->size.x)) * ZZ;
            for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                RR += ((BGRA*)dpos)->R;
                GG += ((BGRA*)dpos)->G;
                BB += ((BGRA*)dpos)->B;
                M = max(((BGRA*)dpos)->R, max(((BGRA*)dpos)->G, max(((BGRA*)dpos)->B, M)));
                dpos += ZZ;
            }
            Progress(vcnf);
        }

        CR = (RR)? (FLOAT)(RR + GG + BB)/(RR*3) : 1.0;
        CG = (GG)? (FLOAT)(RR + GG + BB)/(GG*3) : 1.0;
        CB = (BB)? (FLOAT)(RR + GG + BB)/(BB*3) : 1.0;
        X = (M)? 255.0/(M*max(max(CR, CG), CB)) : 1.0;
        CR *= X;
        CG *= X;
        CB *= X;

        for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
            dpos = src + (vcnf->area.left + (y * curr->size.x)) * ZZ;
            for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                ((BGRA*)dpos)->R *= CR;
                ((BGRA*)dpos)->G *= CG;
                ((BGRA*)dpos)->B *= CB;
                dpos += ZZ;
            }
            Progress(vcnf);
        }
    }
}



inline void MedianFltAddRows(HIST *hist, BYTE *src, RECT *crop, DWORD width) {
    DWORD x;

    src += (crop->left + (crop->bottom * width)) * ZZ;
    for (x = crop->left; x <= crop->right; x++) {
        hist->rpix[((BGRA*)src)->R]++;
        hist->gpix[((BGRA*)src)->G]++;
        hist->bpix[((BGRA*)src)->B]++;
        hist->npix++;
        src += ZZ;
    }
}

inline void MedianFltDelRows(HIST *hist, BYTE *src, RECT *crop, DWORD width) {
    DWORD x;

    src += (crop->left + (crop->top * width)) * ZZ;
    for (x = crop->left; x <= crop->right; x++) {
        hist->rpix[((BGRA*)src)->R]--;
        hist->gpix[((BGRA*)src)->G]--;
        hist->bpix[((BGRA*)src)->B]--;
        hist->npix--;
        src += ZZ;
    }
}

DWORD CALLBACK MedianFlt_thr(TCNF *tcnf) {
    CONF *vcnf = tcnf->conf;
    BYTE *src = curr->prev->bptr;
    BYTE *dst = curr->bptr;
    HIST hist;

    LONG x, y, dpos, mrad = vcnf->data & 0xFF;
    RECT crop;

    for (x = tcnf->dmin; x < tcnf->dmax; x++) {
        crop.left = max(vcnf->area.left, x - mrad);
        crop.right = min(vcnf->area.right - 1, x + mrad);

        hist.npix = 0;
        for (y = 0; y < 256; y++)
            hist.rpix[y] = hist.gpix[y] = hist.bpix[y] = 0;

        for (y = vcnf->area.top; y < vcnf->area.top + mrad; y++)
            if (y < vcnf->area.bottom) {
                crop.bottom = y;
                MedianFltAddRows(&hist, src, &crop, curr->size.x);
            }

        dpos = (x + (vcnf->area.top * curr->size.x)) * ZZ;
        for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
            crop.top = y - mrad;
            crop.bottom = y + mrad;

            if (crop.top >= vcnf->area.top)
                MedianFltDelRows(&hist, src, &crop, curr->size.x);
            if (crop.bottom < vcnf->area.bottom)
                MedianFltAddRows(&hist, src, &crop, curr->size.x);

            DWORD npd2 = (hist.npix + 1) >> 1;
            if (!npd2)
                ((BGRA*)(dst + dpos))->BGRA = ((BGRA*)(src + dpos))->BGRA;
            else {
                DWORD i, s;
                for (i = s = 0; (s += hist.rpix[i]) < npd2; i++);
                ((BGRA*)(dst + dpos))->R = i;
                for (i = s = 0; (s += hist.gpix[i]) < npd2; i++);
                ((BGRA*)(dst + dpos))->G = i;
                for (i = s = 0; (s += hist.bpix[i]) < npd2; i++);
                ((BGRA*)(dst + dpos))->B = i;
            }
            dpos += curr->size.x * ZZ;
        }
        Progress(vcnf);
    }
    return TRUE;
}

void MedianFlt(CONF *vcnf, DWORD iact) {
    if (!iact) {
        HWND hDlg = FindWindowEx(vcnf->hwid, 0, DLG_WNDC, 0);
        vcnf->data = abs((SendMessage(GetDlgItem(hDlg, FUE_GRAD + 1), UDM_GETPOS, 0, 0) & 0xFF) % 101);
    }
    else {
        if (!(vcnf->data & 0xFF)) return;
        vcnf->mprg = 100.0/(vcnf->area.right - vcnf->area.left);
        Parallelize(vcnf, vcnf->area.left, vcnf->area.right, MedianFlt_thr);
    }
}



void AContLvls(CONF *vcnf, DWORD iact) {
    if (!iact) {
        HWND hDlg = FindWindowEx(vcnf->hwid, 0, DLG_WNDC, 0);
        vcnf->data = ((SendMessage(GetDlgItem(hDlg, FCR_ACTR), BM_GETCHECK, 0, 0))? 0x100 : 0)
                |    ((SendMessage(GetDlgItem(hDlg, FCR_ALVL), BM_GETCHECK, 0, 0))? 0x200 : 0)
                | abs((SendMessage(GetDlgItem(hDlg, FCE_CLIP + 1), UDM_GETPOS, 0, 0) & 0xFF) % 51);
    }
    else {
        LONG x, y, dpos, lclr, lclg, lclb, hclr, hclg, hclb;
        BYTE *src = curr->bptr;
        FLOAT coer, coeg, coeb;
        HIST hist = {};

        if (vcnf->data & 0xFF00) {
            hclr = hclg = hclb = 0;
            lclr = lclg = lclb = 255;

            vcnf->mprg = 50.0/(vcnf->area.bottom - vcnf->area.top);

            for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
                dpos = (vcnf->area.left + (y * curr->size.x)) * ZZ;
                for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                    hist.rpix[((BGRA*)(src + dpos))->R]++;
                    hist.gpix[((BGRA*)(src + dpos))->G]++;
                    hist.bpix[((BGRA*)(src + dpos))->B]++;
                    dpos += ZZ;
                }
                Progress(vcnf);
            }

            dpos = (vcnf->area.bottom - vcnf->area.top)*(vcnf->area.right - vcnf->area.left)*(vcnf->data & 0xFF)/100;
            for (x = 0, lclr =  -1; x <= dpos; x += hist.rpix[++lclr]);
            for (x = 0, lclg =  -1; x <= dpos; x += hist.gpix[++lclg]);
            for (x = 0, lclb =  -1; x <= dpos; x += hist.bpix[++lclb]);
            for (x = 0, hclr = 256; x <= dpos; x += hist.rpix[--hclr]);
            for (x = 0, hclg = 256; x <= dpos; x += hist.gpix[--hclg]);
            for (x = 0, hclb = 256; x <= dpos; x += hist.bpix[--hclb]);

            if (vcnf->data & 256) {
                lclr = lclg = lclb = min(lclr, min(lclg, lclb));
                hclr = hclg = hclb = max(hclr, min(hclg, hclb));
            }
            coer = 256.0/(hclr - lclr + 1);
            coeg = 256.0/(hclg - lclg + 1);
            coeb = 256.0/(hclb - lclb + 1);
            for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
                dpos = (vcnf->area.left + (y * curr->size.x)) * ZZ;
                for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                    ((BGRA*)(src + dpos))->R = min(255, (FLOAT)max(0, ((BGRA*)(src + dpos))->R - (LONG)lclr) * coer);
                    ((BGRA*)(src + dpos))->G = min(255, (FLOAT)max(0, ((BGRA*)(src + dpos))->G - (LONG)lclg) * coeg);
                    ((BGRA*)(src + dpos))->B = min(255, (FLOAT)max(0, ((BGRA*)(src + dpos))->B - (LONG)lclb) * coeb);
                    dpos += ZZ;
                }
                Progress(vcnf);
            }
        }
    }
}



DWORD CALLBACK GaussBlurH_thr(TCNF *tcnf) {
    CONF *vcnf = tcnf->conf;
    BYTE *bpl, *bph, *dst = curr->bptr, *src = curr->prev->bptr;
    FLOAT rsum, gsum, bsum, *vect = (FLOAT*)(vcnf->data + sizeof(FLOAT));
    LONG x, y, z, dpos, blur = *((FLOAT*)vcnf->data);

    for (y = tcnf->dmin; y < tcnf->dmax; y++) {
        dpos = (y * curr->size.x) * ZZ;
        for (x = vcnf->area.left; x < vcnf->area.right; x++) {
            rsum = gsum = bsum = 0.0;
            for (z = 1; z <= blur; z++) {
                bpl = src + dpos + ZZ * max(vcnf->area.left,      x - z);
                bph = src + dpos + ZZ * min(vcnf->area.right - 1, x + z);
                rsum += ((((BGRA*)bpl)->R) + (((BGRA*)bph)->R)) * vect[z];
                gsum += ((((BGRA*)bpl)->G) + (((BGRA*)bph)->G)) * vect[z];
                bsum += ((((BGRA*)bpl)->B) + (((BGRA*)bph)->B)) * vect[z];
            }
            bpl = dst + dpos + ZZ * x;
            ((BGRA*)bpl)->R = ((BGRA*)bpl)->R * vect[0] + rsum;
            ((BGRA*)bpl)->G = ((BGRA*)bpl)->G * vect[0] + gsum;
            ((BGRA*)bpl)->B = ((BGRA*)bpl)->B * vect[0] + bsum;
        }
        Progress(vcnf);
    }
    return TRUE;
}

DWORD CALLBACK GaussBlurV_thr(TCNF *tcnf) {
    CONF *vcnf = tcnf->conf;
    BYTE *bpl, *bph, *dst = curr->bptr, *src = curr->next->bptr;
    FLOAT rsum, gsum, bsum, *vect = (FLOAT*)(vcnf->data + sizeof(FLOAT));
    LONG x, y, z, blur = *((FLOAT*)vcnf->data);

    for (x = tcnf->dmin; x < tcnf->dmax; x++) {
        for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
            rsum = gsum = bsum = 0.0;
            for (z = 1; z <= blur; z++) {
                bpl = src + ZZ * (x + curr->size.x * max(vcnf->area.top,        y - z));
                bph = src + ZZ * (x + curr->size.x * min(vcnf->area.bottom - 1, y + z));
                rsum += ((((BGRA*)bpl)->R) + (((BGRA*)bph)->R)) * vect[z];
                gsum += ((((BGRA*)bpl)->G) + (((BGRA*)bph)->G)) * vect[z];
                bsum += ((((BGRA*)bpl)->B) + (((BGRA*)bph)->B)) * vect[z];
            }
            bpl = dst + ZZ * (x + curr->size.x * max(vcnf->area.top, y));
            ((BGRA*)bpl)->R = ((BGRA*)bpl)->R * vect[0] + rsum;
            ((BGRA*)bpl)->G = ((BGRA*)bpl)->G * vect[0] + gsum;
            ((BGRA*)bpl)->B = ((BGRA*)bpl)->B * vect[0] + bsum;
        }
        Progress(vcnf);
    }
    return TRUE;
}

void GaussBlur(CONF *vcnf, DWORD iact) {
    if (!iact) {
        HWND hDlg = FindWindowEx(vcnf->hwid, 0, DLG_WNDC, 0);
        vcnf->data = abs(SendMessage(GetDlgItem(hDlg, FUE_GRAD + 1), UDM_GETPOS, 0, 0) % 101) * 1000;
    }
    else {
        FLOAT *vect, rsum, gsum, bsum;
        LONG x;

        if (!vcnf->data) return;

        gsum = ((FLOAT)vcnf->data / 1000.0) * 3.0;
        bsum = 9.0 / (2.0 * gsum * gsum);
        iact = tr(gsum);
        vect = (FLOAT*)malloc((iact + 2) * sizeof(FLOAT));

        rsum = 0.0;
        vect[1] = 1.0;
        for (x = 1; x <= iact; x++)
            rsum += vect[x + 1] = exp(-x * x * bsum);

        rsum = 1.0 / (rsum + rsum + 1.0);
        for (x = 0; x <= iact; x++)
            vect[x + 1] *= rsum;

        vect[0] = iact;
        vcnf->data = (DWORD)vect;
        vcnf->mprg = 50.0/(vcnf->area.bottom - vcnf->area.top);

        Parallelize(vcnf, vcnf->area.top, vcnf->area.bottom, GaussBlurH_thr);

        DuplicatePict(curr);
        vcnf->cprg = 50.0;
        vcnf->dprg = vcnf->area.right - vcnf->area.left;
        vcnf->mprg = vcnf->cprg/vcnf->dprg;

        Parallelize(vcnf, vcnf->area.left, vcnf->area.right, GaussBlurV_thr);

        FreePictTree(&curr->next);
        free(vect);
    }
}



void SobelEdge(CONF *vcnf, DWORD iact) {
    if (!iact) {
        HWND hDlg = FindWindowEx(vcnf->hwid, 0, DLG_WNDC, 0);
        vcnf->data = ((SendMessage(GetDlgItem(hDlg, FBR_HORZ), BM_GETCHECK, 0, 0))? 1 : 0)
                   | ((SendMessage(GetDlgItem(hDlg, FBR_VERT), BM_GETCHECK, 0, 0))? 2 : 0)
                   | ((SendMessage(GetDlgItem(hDlg, FBC_FREV), BM_GETCHECK, 0, 0))? 4 : 0);
    }
    else {
        LONG x, y, z, frev, dpos;

        PICT *ptmp = MakePict();
        AssignDIB(ptmp, curr->hdib, 2*curr->size.x, -curr->size.y);

        SHORT *dst = (SHORT*)ptmp->bptr;
        BYTE *src = curr->bptr;

        frev = (vcnf->data & 4)? -1 : 1;
        vcnf->mprg = 100.0/((1 + ((vcnf->data & 1)? 2 : 0) + ((vcnf->data & 2)? 2 : 0)) * (vcnf->area.bottom - vcnf->area.top));

        PICT *horz = MakePict();
        AssignDIB(horz, curr->hdib, 2*curr->size.x, -curr->size.y);

        PICT *vert = MakePict();
        AssignDIB(vert, curr->hdib, 2*curr->size.x, -curr->size.y);

        SHORT *hptr = (SHORT*)horz->bptr, *vptr = (SHORT*)vert->bptr;

        if (vcnf->data & 1) {
            for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
                dpos = (vcnf->area.left + 1 + (y * curr->size.x)) * ZZ;
                for (x = vcnf->area.left + 1; x < vcnf->area.right - 1; x++) {
                    dst[dpos + 0] = frev * ((((BGRA*)(src + dpos - ZZ))->B - (LONG)((BGRA*)(src + dpos + ZZ))->B) << 1);
                    dst[dpos + 1] = frev * ((((BGRA*)(src + dpos - ZZ))->G - (LONG)((BGRA*)(src + dpos + ZZ))->G) << 1);
                    dst[dpos + 2] = frev * ((((BGRA*)(src + dpos - ZZ))->R - (LONG)((BGRA*)(src + dpos + ZZ))->R) << 1);
                    dpos += ZZ;
                }
                Progress(vcnf);
            }
            z = curr->size.x * ZZ;
            for (y = vcnf->area.top + 1; y < vcnf->area.bottom - 1; y++) {
                dpos = (vcnf->area.left + (y * curr->size.x)) * ZZ;
                for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                    hptr[dpos + 0] = dst[dpos + 0] + dst[dpos + z + 0] + dst[dpos - z + 0];
                    hptr[dpos + 1] = dst[dpos + 1] + dst[dpos + z + 1] + dst[dpos - z + 1];
                    hptr[dpos + 2] = dst[dpos + 2] + dst[dpos + z + 2] + dst[dpos - z + 2];
                    dpos += ZZ;
                }
                Progress(vcnf);
            }
        }
        if (vcnf->data & 2) {
            for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
                dpos = (vcnf->area.left + 1 + (y * curr->size.x)) * ZZ;
                for (x = vcnf->area.left + 1; x < vcnf->area.right - 1; x++) {
                    dst[dpos + 0] = ((BGRA*)(src + dpos - ZZ))->B + ((BGRA*)(src + dpos + ZZ))->B + (LONG)((((BGRA*)(src + dpos))->B) << 1);
                    dst[dpos + 1] = ((BGRA*)(src + dpos - ZZ))->G + ((BGRA*)(src + dpos + ZZ))->G + (LONG)((((BGRA*)(src + dpos))->G) << 1);
                    dst[dpos + 2] = ((BGRA*)(src + dpos - ZZ))->R + ((BGRA*)(src + dpos + ZZ))->R + (LONG)((((BGRA*)(src + dpos))->R) << 1);
                    dpos += ZZ;
                }
                Progress(vcnf);
            }
            z = curr->size.x * ZZ;
            for (y = vcnf->area.top + 1; y < vcnf->area.bottom - 1; y++) {
                dpos = (vcnf->area.left + (y * curr->size.x)) * ZZ;
                for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                    vptr[dpos + 0] = frev * (dst[dpos + z + 0] - (LONG)dst[dpos - z + 0]);
                    vptr[dpos + 1] = frev * (dst[dpos + z + 1] - (LONG)dst[dpos - z + 1]);
                    vptr[dpos + 2] = frev * (dst[dpos + z + 2] - (LONG)dst[dpos - z + 2]);
                    dpos += ZZ;
                }
                Progress(vcnf);
            }
        }
        if ((vcnf->data & 3) == 1) {
            for (y = vcnf->area.top + 1; y < vcnf->area.bottom - 1; y++) {
                dpos = (vcnf->area.left + (y * curr->size.x)) * ZZ;
                for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                    ((BGRA*)(src + dpos - ZZ))->B = min(255, max(0, hptr[dpos + 0]));
                    ((BGRA*)(src + dpos - ZZ))->G = min(255, max(0, hptr[dpos + 1]));
                    ((BGRA*)(src + dpos - ZZ))->R = min(255, max(0, hptr[dpos + 2]));
                    dpos += ZZ;
                }
                Progress(vcnf);
            }
        }
        else if ((vcnf->data & 3) == 2) {
            for (y = vcnf->area.top + 1; y < vcnf->area.bottom - 1; y++) {
                dpos = (vcnf->area.left + (y * curr->size.x)) * ZZ;
                for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                    ((BGRA*)(src + dpos - ZZ))->B = min(255, max(0, vptr[dpos + 0]));
                    ((BGRA*)(src + dpos - ZZ))->G = min(255, max(0, vptr[dpos + 1]));
                    ((BGRA*)(src + dpos - ZZ))->R = min(255, max(0, vptr[dpos + 2]));
                    dpos += ZZ;
                }
                Progress(vcnf);
            }
        }
        else if ((vcnf->data & 3) == 3) {
            for (y = vcnf->area.top + 1; y < vcnf->area.bottom - 1; y++) {
                dpos = (vcnf->area.left + (y * curr->size.x)) * ZZ;
                for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                    ((BGRA*)(src + dpos - ZZ))->B = min(255, sqrt((LONG)hptr[dpos + 0]*hptr[dpos + 0] + (LONG)vptr[dpos + 0]*vptr[dpos + 0]));
                    ((BGRA*)(src + dpos - ZZ))->G = min(255, sqrt((LONG)hptr[dpos + 1]*hptr[dpos + 1] + (LONG)vptr[dpos + 1]*vptr[dpos + 1]));
                    ((BGRA*)(src + dpos - ZZ))->R = min(255, sqrt((LONG)hptr[dpos + 2]*hptr[dpos + 2] + (LONG)vptr[dpos + 2]*vptr[dpos + 2]));
                    dpos += ZZ;
                }
                Progress(vcnf);
            }
        }
        FreePictTree(&horz);
        FreePictTree(&vert);

        FreePictTree(&ptmp);
    }
}



void CustomLin(CONF *vcnf, DWORD iact) {
    if (!iact) {
        HWND hTxt = GetDlgItem(FindWindowEx(vcnf->hwid, 0, DLG_WNDC, 0), FNE_MATR);
        iact = SendMessage(hTxt, WM_GETTEXTLENGTH, 0, 0);
        vcnf->data = (DWORD)malloc(iact + 2);
        SendMessage(hTxt, WM_GETTEXT, iact + 2, (LPARAM)vcnf->data);
    }
    else {
        LPSTR mstr, mptr;
        LONG dpos, ipos, jpos, i, j, x = 0, y = 0, z = 0, w = 0;
        BYTE *off, *src = curr->prev->bptr, *dst = curr->bptr;
        FLOAT tmpr, tmpg, tmpb, *mtrx;

        if (*(mstr = mptr = (LPSTR)vcnf->data)) {
            do {
                if (*mptr == ',') z++;
                if (*mptr == ';') {
                    y++;
                    x = max(x, z);
                    z = 0;
                }
                while (*mptr == ' ' || *mptr == '\x0A' || *mptr == '\x0D') mptr++;
            } while ((*mstr++ = *mptr++));
            if (mstr[-2] != ';') y++;
            x = max(x, z) + 1;

            mtrx = (FLOAT*)calloc(x * y, sizeof(FLOAT));
            mptr = mstr = (LPSTR)vcnf->data;

            z = 0;
            w = 1;
            while (*mstr) {
                while (*mptr != ',' && *mptr != ';' && *mptr != '\0') mptr++;
                if (*mptr == ';') w = -w;
                if (*mptr)
                    *mptr++ = '\0';
                else
                    mptr = NULL;

                sscanf(mstr, "%f", &mtrx[x * (abs(w) - 1) + z++]);
                if (w < 0 || z >= x) {
                    w = abs(w) + 1;
                    z = 0;
                }
                if (!mptr) break;

                mstr = mptr;
            }

            ipos = (x - 1) >> 1;
            jpos = (y - 1) >> 1;
            vcnf->mprg = 100.0/(vcnf->area.bottom - vcnf->area.top);

            for (w = vcnf->area.top; w < vcnf->area.bottom; w++) {
                dpos = (vcnf->area.left + (w * curr->size.x)) * ZZ;
                for (z = vcnf->area.left; z < vcnf->area.right; z++) {
                    tmpr = tmpg = tmpb = 0.0;
                    for (j = 0; j < y; j++)
                        if ((w + (j - jpos) >= vcnf->area.top) && (w + (j - jpos) < vcnf->area.bottom)) {
                            off = src + dpos + (-ipos + (j - jpos) * curr->size.x) * ZZ;
                            for (i = 0; i < x; i++)
                                if ((z + (i - ipos) >= vcnf->area.left) && (z + (i - ipos) < vcnf->area.right)) {
                                    tmpr += ((BGRA*)(off + i * ZZ))->R * mtrx[x*j + i];
                                    tmpg += ((BGRA*)(off + i * ZZ))->G * mtrx[x*j + i];
                                    tmpb += ((BGRA*)(off + i * ZZ))->B * mtrx[x*j + i];
                                }
                        }
                    ((BGRA*)(dst + dpos))->R = min(255, max(0, tmpr));
                    ((BGRA*)(dst + dpos))->G = min(255, max(0, tmpg));
                    ((BGRA*)(dst + dpos))->B = min(255, max(0, tmpb));
                    dpos += ZZ;
                }
                Progress(vcnf);
            }
            free(mtrx);
        }
        free((HGLOBAL)vcnf->data);
    }
}



DWORD CALLBACK ScalerFltC_thr(TCNF *tcnf) {
    CONF *vcnf = tcnf->conf;
    FLOAT fscl = 1000.0 / (FLOAT)((vcnf->data & 0xFFFF) % 16001);
    LONG x, y, z, w, dpos, spos;
    BYTE *dst, *src;

    dst = curr->bptr;
    src = curr->prev->bptr;
    for (y = tcnf->dmax; y >= tcnf->dmin; y--) {
        w = min(vcnf->area.bottom - 1, (FLOAT)y * fscl + vcnf->area.top);
        dpos = (y * curr->size.x) * ZZ;
        spos = (w * curr->prev->size.x) * ZZ;

        for (x = curr->size.x - 1; x >= 0; x--) {
            z = min(vcnf->area.right - 1, (FLOAT)x * fscl + vcnf->area.left);
            ((BGRA*)(dst + dpos + x * ZZ))->BGRA = ((BGRA*)(src + spos + z * ZZ))->BGRA;
        }
        Progress(vcnf);
    }
    return TRUE;
}

DWORD CALLBACK ScalerFltL_thr(TCNF *tcnf) {
    CONF *vcnf = tcnf->conf;
    FLOAT dirx, diry, revx, revy, fscx, fscy = 1000.0 / (FLOAT)((vcnf->data & 0xFFFF) % 16001);
    LONG x, y, z, w, dpos, spos, sppl;
    BYTE *dst, *src;

    dst = curr->bptr;
    src = curr->prev->bptr;
    fscx = fscy;
    fscx *= (FLOAT)(vcnf->area.bottom - vcnf->area.top - 1) / (FLOAT)(vcnf->area.bottom - vcnf->area.top);
    fscy *= (FLOAT)(vcnf->area.right - vcnf->area.left - 1) / (FLOAT)(vcnf->area.right - vcnf->area.left);

    #define sp00 (FLOAT)((BGRA*)(src + spos + (z + 0) * ZZ))
    #define sp01 (FLOAT)((BGRA*)(src + spos + (z + 1) * ZZ))
    #define sp10 (FLOAT)((BGRA*)(src + sppl + (z + 0) * ZZ))
    #define sp11 (FLOAT)((BGRA*)(src + sppl + (z + 1) * ZZ))
    for (y = tcnf->dmax; y >= tcnf->dmin; y--) {
        diry = (FLOAT)y * fscy;
        w = min(vcnf->area.bottom - 1, diry + vcnf->area.top);
        revy = 1.0 - (diry -= tr(diry));
        dpos = (y * curr->size.x) * ZZ;
        spos = (w * curr->prev->size.x) * ZZ;
        sppl = spos + curr->prev->size.x * ZZ;

        for (x = curr->size.x - 1; x >= 0; x--) {
            dirx = (FLOAT)x * fscx;
            z = min(vcnf->area.right - 1, dirx + vcnf->area.left);
            revx = 1.0 - (dirx -= tr(dirx));
            ((BGRA*)(dst + dpos + x * ZZ))->B = (sp00->B * revx + sp01->B * dirx) * revy + (sp10->B * revx + sp11->B * dirx) * diry;
            ((BGRA*)(dst + dpos + x * ZZ))->G = (sp00->G * revx + sp01->G * dirx) * revy + (sp10->G * revx + sp11->G * dirx) * diry;
            ((BGRA*)(dst + dpos + x * ZZ))->R = (sp00->R * revx + sp01->R * dirx) * revy + (sp10->R * revx + sp11->R * dirx) * diry;
        }
        Progress(vcnf);
    }
    #undef sp11
    #undef sp10
    #undef sp01
    #undef sp00
    return TRUE;
}

void ScalerFlt(CONF *vcnf, DWORD iact) {
    if (!iact) {
        HWND hDlg = FindWindowEx(vcnf->hwid, 0, DLG_WNDC, 0);
        vcnf->data = ((SendMessage(GetDlgItem(hDlg, FSR_BLIN), BM_GETCHECK, 0, 0))? 0x20000 : 0)
                |    ((SendMessage(GetDlgItem(hDlg, FSR_NONE), BM_GETCHECK, 0, 0))? 0x10000 : 0)
                | abs((SendMessage(GetDlgItem(hDlg, FSE_SCLP + 1), UDM_GETPOS, 0, 0) % 1601) * 10);
    }
    else {
        FLOAT fscl = (FLOAT)((vcnf->data & 0xFFFF) % 16001) / 1000.0;
        LONG x, y;

        if (fscl == 0.0) return;
        if (fscl == 1.0) {
            CutPict(vcnf);
            return;
        }
        curr = curr->prev;
        FreePictTree(&curr->next);
        curr->next = MakePict();
        curr->next->prev = curr;
        curr = curr->next;

        x =  fscl * (vcnf->area.right - vcnf->area.left);
        y = -fscl * (vcnf->area.bottom - vcnf->area.top);
        if (x && y)
            AssignDIB(curr, curr->prev->hdib, x, y);
        else {
            AssignDIB(curr, curr->prev->hdib, 1, -1);
            return;
        }
        vcnf->mprg = 100.0/curr->size.y;
        if (vcnf->data & 0x10000)
            Parallelize(vcnf, 0, curr->size.y - 1, ScalerFltC_thr);
        else if (vcnf->data & 0x20000)
            Parallelize(vcnf, 0, curr->size.y - 1, ScalerFltL_thr);
    }
}



DWORD CALLBACK RotateFlt_thr(TCNF *tcnf) {
    CONF *vcnf = tcnf->conf;
    FLOAT xtmp, ytmp, xvec, yvec, xcnt, ycnt, xadd, yadd,
          qsin, qcos, co00, co01, co10, co11;
    LONG x, y, dpos, qwid, qhei, xsrc, ysrc;
    BYTE *dst, *src;

    dpos = (vcnf->data & 0xFF);
    if (!(vcnf->data & 0x100)) dpos = -dpos;
    qcos = (FLOAT)dpos * M_PI/180.0;
    qsin = sin(qcos);
    qcos = cos(qcos);
    x = vcnf->area.right - vcnf->area.left;
    y = vcnf->area.bottom - vcnf->area.top;
    qwid = tr(x / 2.0);
    qhei = tr(y / 2.0);

    #define sr00 (FLOAT)((BGRA*)(src + ((xsrc + 0) + (ysrc + 0) * curr->prev->size.x) * ZZ))
    #define sr01 (FLOAT)((BGRA*)(src + ((xsrc + 0) + (ysrc + 1) * curr->prev->size.x) * ZZ))
    #define sr10 (FLOAT)((BGRA*)(src + ((xsrc + 1) + (ysrc + 0) * curr->prev->size.x) * ZZ))
    #define sr11 (FLOAT)((BGRA*)(src + ((xsrc + 1) + (ysrc + 1) * curr->prev->size.x) * ZZ))

    src = curr->prev->bptr;
    dst = curr->bptr + (tcnf->dmax + 1) * curr->size.x * ZZ;
    xcnt = (tr(curr->size.x / 2.0) - 0.5) * qcos - (tr(curr->size.y / 2.0) + 0.5) * qsin + 0.5;
    ycnt = (tr(curr->size.x / 2.0) - 0.5) * qsin + (tr(curr->size.y / 2.0) - 0.5) * qcos + 0.5;

    for (y = tcnf->dmax; y >= tcnf->dmin; y--) {
        xadd = - (FLOAT)y * qsin - xcnt;
        yadd = + (FLOAT)y * qcos - ycnt;

        for (x = curr->size.x - 1; x >= 0; x--) {
            dst -= ZZ;

            xvec = (FLOAT)x * qcos + xadd;
            yvec = (FLOAT)x * qsin + yadd;
            if (xvec >= 0.0) xvec += 1.0;
            if (yvec >= 0.0) yvec += 1.0;
            xsrc = qwid + vcnf->area.left + tr(xvec);
            ysrc = qhei + vcnf->area.top  + tr(yvec);

            if (abs(tr(xvec)) >= qwid || abs(tr(yvec)) >= qhei)
                ((BGRA*)dst)->BGRA = CLR_BACK;
            else {
                xvec += qwid;
                yvec += qhei;
                xtmp = xvec - tr(xvec);
                ytmp = yvec - tr(yvec);
                xvec = 1.0 - xtmp;
                yvec = 1.0 - ytmp;

                co00 = xvec * yvec;
                co10 = xtmp * yvec;
                co01 = xvec * ytmp;
                co11 = xtmp * ytmp;

                if (xsrc < vcnf->area.right - 1 && ysrc < vcnf->area.bottom - 1) {
                    ((BGRA*)dst)->R = sr00->R * co00 + sr10->R * co10 + sr01->R * co01 + sr11->R * co11;
                    ((BGRA*)dst)->G = sr00->G * co00 + sr10->G * co10 + sr01->G * co01 + sr11->G * co11;
                    ((BGRA*)dst)->B = sr00->B * co00 + sr10->B * co10 + sr01->B * co01 + sr11->B * co11;
                }
                else if (xsrc < vcnf->area.right - 1) {
                    ((BGRA*)dst)->R = sr00->R * co00 + sr10->R * co10 + ((BGRA)(DWORD)CLR_BACK).R * ytmp;
                    ((BGRA*)dst)->G = sr00->G * co00 + sr10->G * co10 + ((BGRA)(DWORD)CLR_BACK).G * ytmp;
                    ((BGRA*)dst)->B = sr00->B * co00 + sr10->B * co10 + ((BGRA)(DWORD)CLR_BACK).B * ytmp;
                }
                else if (ysrc < vcnf->area.bottom - 1) {
                    ((BGRA*)dst)->R = sr00->R * co00 + sr01->R * co01 + ((BGRA)(DWORD)CLR_BACK).R * xtmp;
                    ((BGRA*)dst)->G = sr00->G * co00 + sr01->G * co01 + ((BGRA)(DWORD)CLR_BACK).G * xtmp;
                    ((BGRA*)dst)->B = sr00->B * co00 + sr01->B * co01 + ((BGRA)(DWORD)CLR_BACK).B * xtmp;
                }
                else {
                    co10 *= ytmp;
                    ((BGRA*)dst)->R = sr00->R * co00 + ((BGRA)(DWORD)CLR_BACK).R * co10;
                    ((BGRA*)dst)->G = sr00->G * co00 + ((BGRA)(DWORD)CLR_BACK).G * co10;
                    ((BGRA*)dst)->B = sr00->B * co00 + ((BGRA)(DWORD)CLR_BACK).B * co10;
                }
            }
        }
        Progress(vcnf);
    }

    #undef sr11
    #undef sr10
    #undef sr01
    #undef sr00

    return TRUE;
}

void RotateFlt(CONF *vcnf, DWORD iact) {
    if (!iact) {
        HWND hDlg = FindWindowEx(vcnf->hwid, 0, DLG_WNDC, 0);
        SHORT angl = SendMessage(GetDlgItem(hDlg, FTE_ANGL + 1), UDM_GETPOS, 0, 0);
        vcnf->data = ((SendMessage(GetDlgItem(hDlg, FTC_SIZE), BM_GETCHECK, 0, 0))? 0x200 : 0)
                   | (abs(angl) % 181) | ((angl < 0)? 0x100 : 0);
    }
    else {
        LONG x, y, dpos;
        FLOAT qsin, qcos;
        BYTE *dst, *src;

        dpos = (vcnf->data & 0xFF);
        if (!dpos) {
            CutPict(vcnf);
            return;
        }

        curr = curr->prev;
        FreePictTree(&curr->next);
        curr->next = MakePict();
        curr->next->prev = curr;
        curr = curr->next;

        if (dpos == 180 || vcnf->data & 0x200) {
            curr->size.x = vcnf->area.right - vcnf->area.left;
            curr->size.y = vcnf->area.bottom - vcnf->area.top;
        }
        else if (dpos == 90) {
            curr->size.x = vcnf->area.bottom - vcnf->area.top;
            curr->size.y = vcnf->area.right - vcnf->area.left;
        }

        if (!(vcnf->data & 0x100)) dpos = -dpos;
        qcos = (FLOAT)dpos * M_PI/180.0;
        qsin = sin(qcos);
        qcos = cos(qcos);

        x = vcnf->area.right - vcnf->area.left;
        y = vcnf->area.bottom - vcnf->area.top;
        if (!curr->size.x) {
            curr->size.x = tr((FLOAT)x * fabs(qcos) + (FLOAT)y * fabs(qsin)) + 1;
            curr->size.y = tr((FLOAT)x * fabs(qsin) + (FLOAT)y * fabs(qcos)) + 1;
            curr->size.x &= -2;
            curr->size.y &= -2;
        }

        AssignDIB(curr, curr->prev->hdib, curr->size.x, -curr->size.y);
        src = curr->prev->bptr;
        dst = curr->bptr + curr->size.y * curr->size.x * ZZ;
        vcnf->mprg = 100.0/curr->size.y;

        if (abs(dpos) == 180) {
            src += (vcnf->area.left + (vcnf->area.top * curr->prev->size.x)) * ZZ;
            for (y = curr->size.y; y > 0; y--) {
                for (x = 0; x < curr->size.x; x++) {
                    dst -= ZZ;
                    ((BGRA*)dst)->BGRA = ((BGRA*)(src + x * ZZ))->BGRA;
                }
                src += curr->prev->size.x * ZZ;
                Progress(vcnf);
            }
            return;
        }
        if (!(vcnf->data & 0x200)) {
            if (dpos == 90) {
                for (y = vcnf->area.left; y < vcnf->area.right; y++) {
                    for (x = vcnf->area.bottom - 1; x >= vcnf->area.top; x--) {
                        dst -= ZZ;
                        ((BGRA*)dst)->BGRA = ((BGRA*)(src + (y + x * curr->prev->size.x) * ZZ))->BGRA;
                    }
                    Progress(vcnf);
                }
                return;
            }
            if (dpos == -90) {
                for (y = vcnf->area.right - 1; y >= vcnf->area.left; y--) {
                    for (x = vcnf->area.top; x < vcnf->area.bottom; x++) {
                        dst -= ZZ;
                        ((BGRA*)dst)->BGRA = ((BGRA*)(src + (y + x * curr->prev->size.x) * ZZ))->BGRA;
                    }
                    Progress(vcnf);
                }
                return;
            }
        }

        Parallelize(vcnf, 0, curr->size.y - 1, RotateFlt_thr);
    }
}



void BrkGlsFlt(CONF *vcnf, DWORD iact) {
    if (!iact) {
        HWND hDlg = FindWindowEx(vcnf->hwid, 0, DLG_WNDC, 0);
        vcnf->data = abs(SendMessage(GetDlgItem(hDlg, FUE_GRAD + 1), UDM_GETPOS, 0, 0) % 101);
    }
    else {
        BYTE *dst = curr->bptr, *src = curr->prev->bptr;
        LONG x, y, xrnd, yrnd, diam, offs, dpos;

        offs = -vcnf->data;
        diam = 1 + (vcnf->data << 1);
        vcnf->mprg = 100.0/(vcnf->area.bottom - vcnf->area.top);
        for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
            dpos = (vcnf->area.left + (y * curr->size.x)) * ZZ;
            for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                xrnd = (rand() % diam) + offs + x;
                yrnd = (rand() % diam) + offs + y;
                xrnd = min(vcnf->area.right - 1, max(vcnf->area.left, xrnd));
                yrnd = min(vcnf->area.bottom - 1, max(vcnf->area.top, yrnd));
                ((BGRA*)(dst + dpos))->BGRA = ((BGRA*)(src + (xrnd + (yrnd * curr->size.x)) * ZZ))->BGRA;
                dpos += ZZ;
            }
            Progress(vcnf);
        }
    }
}



void SineWaves(CONF *vcnf, DWORD iact) {
    if (!iact) {
        HWND hDlg = FindWindowEx(vcnf->hwid, 0, DLG_WNDC, 0);
        vcnf->data = ((SendMessage(GetDlgItem(hDlg, FWR_HORZ), BM_GETCHECK, 0, 0))? 0x100 : 0)
                |    ((SendMessage(GetDlgItem(hDlg, FWR_VERT), BM_GETCHECK, 0, 0))? 0x200 : 0)
                | abs((SendMessage(GetDlgItem(hDlg, FWE_WLEN + 1), UDM_GETPOS, 0, 0) & 0xFF) % 201);
    }
    else {
        BYTE *dst = curr->bptr, *src = curr->prev->bptr;
        LONG x, y, offl, offr, lmbh = (vcnf->data & 0xFF) << 1;
        FLOAT coer, coel, farg = M_PI / lmbh;

        if (lmbh < 4) return;

        if (vcnf->data & 0x100) {
            src += vcnf->area.left * ZZ;
            dst += vcnf->area.left * ZZ;
            vcnf->mprg = 100.0/(vcnf->area.right - vcnf->area.left);

            for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                coel = (FLOAT)(lmbh >> 1) * sin((x - vcnf->area.left + 0.5) * farg);
                coer = coel - tr(coel);

                if (coer < 0.0)
                    coer += 1.0;

                if (coer < 0.5) {
                    coel -= 1.0;
                    coer += 1.0;
                }
                else if (coer == 0.5)
                    coer = 2.0;

                coer /= 2.0;
                offl = tr(coel) % (vcnf->area.bottom - vcnf->area.top);
                coel = 1.0 - coer;
                if (offl < 0)
                    offl += (vcnf->area.bottom - vcnf->area.top);

                offr = (offl + 1) % (vcnf->area.bottom - vcnf->area.top);

                for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
                    if (y + offl < vcnf->area.top)
                        offl += (vcnf->area.bottom - vcnf->area.top);
                    if (y + offl >= vcnf->area.bottom)
                        offl -= (vcnf->area.bottom - vcnf->area.top);

                    offr = (offl + 1) % (vcnf->area.bottom - vcnf->area.top);

                    if (y + offr < vcnf->area.top)
                        offr += (vcnf->area.bottom - vcnf->area.top);
                    if (y + offr >= vcnf->area.bottom)
                        offr -= (vcnf->area.bottom - vcnf->area.top);

                    ((BGRA*)(dst + y * curr->size.x * ZZ))->B = (FLOAT)((BGRA*)(src + (y + offl) * curr->size.x * ZZ))->B * coel + (FLOAT)((BGRA*)(src + (y + offr) * curr->size.x * ZZ))->B * coer;
                    ((BGRA*)(dst + y * curr->size.x * ZZ))->G = (FLOAT)((BGRA*)(src + (y + offl) * curr->size.x * ZZ))->G * coel + (FLOAT)((BGRA*)(src + (y + offr) * curr->size.x * ZZ))->G * coer;
                    ((BGRA*)(dst + y * curr->size.x * ZZ))->R = (FLOAT)((BGRA*)(src + (y + offl) * curr->size.x * ZZ))->R * coel + (FLOAT)((BGRA*)(src + (y + offr) * curr->size.x * ZZ))->R * coer;
                }
                src += ZZ;
                dst += ZZ;
                Progress(vcnf);
            }
        }
        else if (vcnf->data & 0x200) {
            src += vcnf->area.top * curr->size.x * ZZ;
            dst += vcnf->area.top * curr->size.x * ZZ;
            vcnf->mprg = 100.0/(vcnf->area.bottom - vcnf->area.top);

            for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
                coel = (FLOAT)(lmbh >> 1) * sin((y - vcnf->area.top + 0.5) * farg);
                coer = coel - tr(coel);

                if (coer < 0.0)
                    coer += 1.0;

                if (coer < 0.5) {
                    coel -= 1.0;
                    coer += 1.0;
                }
                else if (coer == 0.5)
                    coer = 2.0;

                coer /= 2.0;
                offl = tr(coel) % (vcnf->area.right - vcnf->area.left);
                coel = 1.0 - coer;
                if (offl < 0)
                    offl += (vcnf->area.right - vcnf->area.left);

                offr = (offl + 1) % (vcnf->area.right - vcnf->area.left);

                for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                    if (x + offl < vcnf->area.left)
                        offl += (vcnf->area.right - vcnf->area.left);
                    if (x + offl >= vcnf->area.right)
                        offl -= (vcnf->area.right - vcnf->area.left);

                    offr = (offl + 1) % (vcnf->area.right - vcnf->area.left);

                    if (x + offr < vcnf->area.left)
                        offr += (vcnf->area.right - vcnf->area.left);
                    if (x + offr >= vcnf->area.right)
                        offr -= (vcnf->area.right - vcnf->area.left);

                    ((BGRA*)(dst + x * ZZ))->B = (FLOAT)((BGRA*)(src + (x + offl) * ZZ))->B * coel + (FLOAT)((BGRA*)(src + (x + offr) * ZZ))->B * coer;
                    ((BGRA*)(dst + x * ZZ))->G = (FLOAT)((BGRA*)(src + (x + offl) * ZZ))->G * coel + (FLOAT)((BGRA*)(src + (x + offr) * ZZ))->G * coer;
                    ((BGRA*)(dst + x * ZZ))->R = (FLOAT)((BGRA*)(src + (x + offl) * ZZ))->R * coel + (FLOAT)((BGRA*)(src + (x + offr) * ZZ))->R * coer;
                }
                src += curr->size.x * ZZ;
                dst += curr->size.x * ZZ;
                Progress(vcnf);
            }
        }
    }
}



void FadToGray(CONF *vcnf, DWORD iact) {
    if (!iact) {
        HWND hDlg = FindWindowEx(vcnf->hwid, 0, DLG_WNDC, 0);
        vcnf->data = ((SendMessage(GetDlgItem(hDlg, FFR_LUMI), BM_GETCHECK, 0, 0))? 1 : 0)
                   | ((SendMessage(GetDlgItem(hDlg, FFR_MEAN), BM_GETCHECK, 0, 0))? 2 : 0);
    }
    else {
        BYTE *dpos, *src = curr->bptr;
        BGRA c = {};
        LONG x, y;

        #define pix ((BGRA*)(dpos + x * ZZ))
        vcnf->mprg = 100.0/(vcnf->area.bottom - vcnf->area.top);
        if (vcnf->data & 1) {
            for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
                dpos = src + y * curr->size.x * ZZ;
                for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                    c.R = c.G = c.B = (FLOAT)pix->R * 0.2126 + (FLOAT)pix->G * 0.7152 + (FLOAT)pix->B * 0.0722;
                    ((BGRA*)(dpos + x * ZZ))->BGRA = c.BGRA;
                }
                Progress(vcnf);
            }
        }
        else if (vcnf->data & 2) {
            for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
                dpos = src + y * curr->size.x * ZZ;
                for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                    c.R = c.G = c.B = (pix->R + pix->G + pix->B) / 3;
                    ((BGRA*)(dpos + x * ZZ))->BGRA = c.BGRA;
                }
                Progress(vcnf);
            }
        }
        #undef pix
    }
}



void GlassTile(CONF *vcnf, DWORD iact) {
    if (!iact) {
        HWND hDlg = FindWindowEx(vcnf->hwid, 0, DLG_WNDC, 0);
        vcnf->data = abs(SendMessage(GetDlgItem(hDlg, FIE_NTIL + 1), UDM_GETPOS, 0, 0) % 101);
    }
    else {
        BYTE *src = curr->prev->bptr, *dst = curr->bptr;
        LONG spos, dpos, qwid, qhei,
             xoff, yoff, xmid, ymid,
             x, y, z, w;

        if (vcnf->data < 2) return;

        vcnf->mprg = 100.0/(vcnf->area.bottom - vcnf->area.top);
        qwid = (vcnf->area.right - vcnf->area.left)/vcnf->data;
        qhei = (vcnf->area.bottom - vcnf->area.top)/vcnf->data;
        if (qwid & 1) qwid++;
        if (qhei & 1) qhei++;

        if (!qwid || !qhei) return;

        yoff = ymid = 0;
        for (y = vcnf->area.top; y < vcnf->area.bottom; y++) {
            w = max(vcnf->area.top, min(vcnf->area.bottom - 1, vcnf->area.top + ymid + yoff * 2));

            spos = w * curr->size.x * ZZ;
            dpos = y * curr->size.x * ZZ;

            xoff = xmid = 0;
            for (x = vcnf->area.left; x < vcnf->area.right; x++) {
                z = max(vcnf->area.left, min(vcnf->area.right - 1, vcnf->area.left + xmid + xoff * 2));
                ((BGRA*)(dst + dpos + x * ZZ))->BGRA = ((BGRA*)(src + spos + z * ZZ))->BGRA;

                if (++xoff << 1 == qwid) {
                    xmid += qwid;
                    xoff -= qwid;
                }
            }
            if (++yoff << 1 == qhei) {
                ymid += qhei;
                yoff -= qhei;
            }
            Progress(vcnf);
        }
    }
}



LPSTR greeting = "Выберите фильтр из списка...";
FLTR allFilters[] = {
    {DLG_FROT, FTC_SIZE, -180,  180, RotateFlt, "Поворот"},
    {DLG_FSCL, FSR_BLIN,    1, 1600, ScalerFlt, "Масштабирование"},
    {DLG_FCTR, FCR_ACTR,    0,   50, AContLvls, "Автоконтраст/автоуровни"},
    {DLG_FSBL, FBR_VERT,    1, 1600, SobelEdge, "Нахождение краёв (Собель)"},
    {DLG_FGAU,        0,    1,  100, GaussBlur, "Размывка по Гауссу"},
    {DLG_FLIN,        0,    0,    0, CustomLin, "Линейное преобразование"},
    {DLG_FGAU,        0,    1,  100, MedianFlt, "Медианный фильтр"},
    {DLG_FFTG, FFR_LUMI,    0,    0, FadToGray, "Оттенки серого"},
    {DLG_FEMP,        0,    0,    0, GrayWorld, "\"Серый мир\""},
    {DLG_FWAV, FWR_HORZ,    2,  200, SineWaves, "\"Волны\""},
    {DLG_FGAU,        0,    1,  100, BrkGlsFlt, "\"Стекло\""},
    {DLG_FTIL,        0,    2,  100, GlassTile, "\"Плитка\""}
};



void MakeConf(CONF **root, HWND hwid, FLTRFUNC *func) {
    CONF *temp = *root;
    *root = (CONF*)calloc(1, sizeof(CONF));
    (*root)->func = func;
    (*root)->next = temp;
    (*root)->hwid = hwid;
    (*root)->mark = CreateSolidBrush((COLORREF)(((rand()*rand()) % 0x1000000) | 0x808080));
}



CONF *FindConf(HWND hwid) {
    CONF *temp = conf;
    while (temp) {
        if (temp->hwid == hwid)
            break;
        temp = temp->next;
    }
    return temp;
}



void FreeConf(HWND hwid) {
    CONF *prev = NULL, *temp = conf;
    while (temp) {
        if (temp->hwid == hwid) {
            DeleteObject(temp->mark);
            if (!prev)
                conf = temp->next;
            else
                prev->next = temp->next;
            free(temp);
            return;
        }
        prev = temp;
        temp = temp->next;
    }
}



PICT *MakePict() {
    PICT *retn = (PICT*)calloc(1, sizeof(PICT));
    retn->devc = CreateCompatibleDC(NULL);
    return retn;
}



void AssignDIB(PICT *pict, HBITMAP hbm, LONG w, LONG h) {
    BITMAPINFO bmi = {};

    if ((w && h) || (!w && !h)) {
        bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
        GetDIBits(pict->devc, hbm, 0, 0, NULL, &bmi, DIB_RGB_COLORS);
        bmi.bmiHeader.biBitCount = 8*ZZ;
        bmi.bmiHeader.biCompression = BI_RGB;

        if (w || h) {
            bmi.bmiHeader.biWidth  = (bmi.bmiHeader.biWidth  > 0)? abs(w) : -abs(w);
            bmi.bmiHeader.biHeight = (bmi.bmiHeader.biHeight > 0)? abs(h) : -abs(h);
        }
        pict->size.x = abs(bmi.bmiHeader.biWidth);
        pict->size.y = abs(bmi.bmiHeader.biHeight);

        pict->hdib = CreateDIBSection(pict->devc, &bmi, DIB_RGB_COLORS, &pict->bptr, 0, 0);
        if (h >= 0)
            GetDIBits(pict->devc, hbm, 0, pict->size.y, pict->bptr, &bmi, DIB_RGB_COLORS);

        SelectObject(pict->devc, pict->hdib);
    }
    else {
        pict->size.x = pict->size.y = 0;
        pict->hdib = NULL;
    }
}



void FreePictTree(PICT **root) {
    while (*root) {
        PICT *next = (*root)->next;
        DeleteDC((*root)->devc);
        DeleteObject((*root)->hdib);
        free(*root);
        *root = next;
    }
}



void DuplicatePict(PICT *pict) {
    FreePictTree(&pict->next);
    pict->next = MakePict();
    AssignDIB(pict->next, pict->hdib, 0, 0);
    pict->next->prev = pict;
    pict->next->size = pict->size;
}



BOOL LoadRootPict(LPSTR fnm, HWND hDlg) {
    BOOL ret = TRUE;
    IPicture *pic;
    IStream *pst;
    HGLOBAL ptr;
    HANDLE hfl;
    DWORD fsz;
    HWND hwc;

    if (hDlg) {
        if (IsWindowEnabled(hwc = GetDlgItem(hDlg, MCB_FLTR))) {
            SendMessage(hwc, CB_INSERTSTRING, 0, (LPARAM)greeting);
            SendMessage(hwc, CB_SETCURSEL, 0, 0);
            EnableWindow(hwc, FALSE);
        }
        EnableWindow(GetDlgItem(hDlg, MBT_DONE), FALSE);
        EnableWindow(GetDlgItem(hDlg, MBT_SAVE), FALSE);
        EnableWindow(fwdButton, FALSE);
        EnableWindow(bwdButton, FALSE);
        flgPaint = FALSE;
    }
    hfl = CreateFile(fnm, GENERIC_READ, 0, 0, OPEN_EXISTING, 0, 0);
    fsz = GetFileSize(hfl, NULL);
    ptr = NULL;

    if ((fsz != (DWORD)-1) &&
        (hfl != INVALID_HANDLE_VALUE) &&
        (ptr  = GlobalAlloc(GMEM_FIXED, fsz))) {

        if (SUCCEEDED(CreateStreamOnHGlobal(ptr, TRUE, &pst)) &&
            ReadFile(hfl, GlobalLock(ptr), fsz, (LPDWORD)&hwc, 0)) {

            FreePictTree(&root);
            root = MakePict();
            OleLoadPicture(pst, fsz, FALSE, (IID*)&IID_IPicture, (void**)&pic);
            curr = root;

            if (pic) {
                pic->lpVtbl->get_Handle(pic, (OLE_HANDLE*)&root->hdib);
                if (root->hdib) {
                    AssignDIB(root, root->hdib, 0, 0);
                    pic->lpVtbl->Release(pic);

                    cpos.x = cpos.y = 0;
                    if (hDlg) {
                        flgPaint = TRUE;
                        InvalidateRect(mainSurface, NULL, FALSE);
                        EnableWindow(GetDlgItem(hDlg, MBT_DONE), TRUE);
                        EnableWindow(GetDlgItem(hDlg, MBT_SAVE), TRUE);
                        EnableWindow(hwc = GetDlgItem(hDlg, MCB_FLTR), TRUE);
                        SendMessage(hwc, CB_DELETESTRING, 0, 0);
                        SendMessage(hwc, CB_SETCURSEL, 0, 0);
                    }
                }
            }
            if (!root->hdib) {
                FreePictTree(&root);
                ret = FALSE;
                if (hDlg) {
                    InvalidateRect(mainSurface, NULL, FALSE);
                    MessageBoxIcon(hDlg, "Неподдерживаемый формат!", NULL, MB_OK | MB_USERICON, (LPSTR)ICN_MAIN);
                }
            }
        }
        else {
            ret = FALSE;
            if (hDlg)
                MessageBoxIcon(hDlg, "Ошибка библиотеки OLE!", NULL, MB_OK | MB_USERICON, (LPSTR)ICN_MAIN);
        }
        if (pst) pst->lpVtbl->Release(pst);
        GlobalUnlock(ptr);
    }
    else {
        ret = FALSE;
        if (hDlg)
            MessageBoxIcon(hDlg, "Ошибка чтения файла!", NULL, MB_OK | MB_USERICON, (LPSTR)ICN_MAIN);
    }
    GlobalFree(ptr);
    CloseHandle(hfl);
    return ret;
}



BOOL SaveCurrPict(LPSTR fnm, HWND hDlg) {
    DWORD ret, flg = flgPaint;
    LONG fsz, wsz;
    BMPH bmh = {};
    BYTE *ptr;

    flgPaint = FALSE;
    bmh.pln = 1;
    bmh.img = 0x28;
    bmh.bmp = 0x4D42;
    bmh.bmw = curr->size.x;
    bmh.bmh = curr->size.y;
    bmh.off = sizeof(BMPH);
    bmh.bpp = ZZ * 8;
    bmh.bsz = ZZ * bmh.bmw * bmh.bmh + bmh.off;

    HANDLE fbm = CreateFile(fnm, GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    WriteFile(fbm, &bmh, bmh.off, &ret, NULL);

    if (!ret) {
        if (hDlg)
            MessageBoxIcon(hDlg, "Ошибка записи заголовка!", NULL, MB_OK | MB_USERICON, (LPSTR)ICN_MAIN);
    }
    else {
        fsz = wsz = bmh.bsz - bmh.off;
        ptr = curr->bptr;

        while (fsz > 0) {
            WriteFile(fbm, ptr, min(fsz, wsz), &ret, NULL);
            if (!ret) {
                wsz = (wsz >> 1) + 1;
                if (wsz < sizeof(BMPH)) {
                    ret = fsz = 0;
                    if (hDlg)
                        MessageBoxIcon(hDlg, "Всё перепробовал, - не записывается! =(", NULL, MB_OK | MB_USERICON, (LPSTR)ICN_MAIN);
                }
            }
            else {
                fsz -= ret;
                ptr +=ret;
            }
        }
    }
    CloseHandle(fbm);
    flgPaint = flg;
    return ret;
}



void NotifyBusy(HWND hDlg) {
    MessageBoxIcon(hDlg, "Пожалуйста подождите, один или несколько\nфильтров в данный момент активны!", NULL, MB_OK | MB_USERICON, (LPSTR)ICN_MAIN);
}



DWORD CALLBACK FilterThread(LPVOID lpPrm) {
    ((CONF*)lpPrm)->func((CONF*)lpPrm, TRUE);
    SendMessage(progBar, PBM_SETPOS, 0, 0);
    flgPaint = TRUE;
    EnableWindow(bwdButton, TRUE);
    SetFocus(bwdButton);
    SendMessage(((CONF*)lpPrm)->hwid, WM_CLOSE, 0, 0);
    InvalidateRect(mainSurface, NULL, FALSE);
    return TRUE;
}



BOOL CALLBACK FltrProc(HWND hDlg, UINT uMsg, WPARAM wPrm, LPARAM lPrm) {
    #define CLR_MGRN 5
    switch (uMsg) {
        case WM_INITDIALOG: {
            RECT rcw;
            GetClientRect(hDlg, &rcw);
            rcw.left = rcw.right - ((rcw.right - rcw.left) >> 3);
            CreateWindowEx(0, WC_STATIC, NULL, SS_SUNKEN | WS_VISIBLE | WS_CHILD, rcw.left, rcw.top + CLR_MGRN, rcw.right - rcw.left - CLR_MGRN, rcw.bottom - rcw.top - 2*CLR_MGRN, hDlg, (HMENU)DLG_MAIN, base, NULL);
            return TRUE;
        }

        case WM_CTLCOLORSTATIC:
            if ((HWND)lPrm == GetDlgItem(hDlg, DLG_MAIN))
                return (DWORD)FindConf(GetParent(hDlg))->mark;
            return FALSE;

        case WM_CLOSE:
            DestroyWindow(GetDlgItem(hDlg, DLG_MAIN));
            EndDialog(hDlg, 0);

        default:
            return FALSE;
    }
    #undef CLR_MRGN
}



BOOL CALLBACK ConfProc(HWND hDlg, UINT uMsg, WPARAM wPrm, LPARAM lPrm) {
    switch (uMsg) {
        case WM_INITDIALOG: {
            HWND hCld, hUpd;

            MakeConf(&conf, hDlg, allFilters[lPrm].func);
            SendMessage(hDlg, WM_SETTEXT, 0, (LPARAM)allFilters[lPrm].name);
            SendMessage(hDlg, WM_SETICON, ICON_BIG, (LPARAM)LoadIcon(base, MAKEINTRESOURCE(ICN_MAIN)));
            hUpd = CreateDialogParam(base, (LPSTR)(DWORD)allFilters[lPrm].tdlg, hDlg, (DLGPROC)FltrProc, 0);
            if (allFilters[lPrm].bdef)
                SendMessage(GetDlgItem(hUpd, allFilters[lPrm].bdef), BM_SETCHECK, BST_CHECKED, 0);
            if (allFilters[lPrm].liml != allFilters[lPrm].limh) {
                hCld = FindWindowEx(hUpd, 0, WC_EDIT, 0);
                hUpd = CreateWindowEx(0, UPDOWN_CLASS, NULL, UDS_HORZ | UDS_ALIGNRIGHT | UDS_ARROWKEYS | UDS_WRAP | UDS_HOTTRACK | UDS_SETBUDDYINT | UDS_NOTHOUSANDS | WS_CHILD | WS_VISIBLE | WS_BORDER, 0, 0, 0, 0, hUpd, (HMENU)(GetDlgCtrlID(hCld) + 1), base, NULL);
                SendMessage(hUpd, UDM_SETBUDDY, (WPARAM)hCld, 0);
                SendMessage(hUpd, UDM_SETRANGE, 0, MAKELPARAM(allFilters[lPrm].limh, allFilters[lPrm].liml));
                SendMessage(hUpd, UDM_SETPOS, 0, 0);
            }
            return TRUE;
        }

        case WM_COMMAND:
            switch (LOWORD(wPrm)) {
                case TCB_RECT:
                    if (!idRect && flgPaint) {
                        SendMessage((HWND)lPrm, BM_SETCHECK, BST_CHECKED, 0);
                        idRect = hDlg;
                    }
                    return FALSE;

                case TBT_FIRE: {
                    if (flgPaint) {
                        flgPaint = FALSE;

                        CONF *temp = FindConf(hDlg);
                        RECT real = temp->area;

                        EnableWindow(fwdButton, FALSE);
                        DuplicatePict(curr);
                        curr = curr->next;

                        if (real.left > real.right) {
                            temp->area.left  = real.right;
                            temp->area.right = real.left;
                        }
                        if (real.top > real.bottom) {
                            temp->area.top    = curr->size.y - real.top;
                            temp->area.bottom = curr->size.y - real.bottom;
                        }
                        else if (real.top < real.bottom) {
                            temp->area.top    = curr->size.y - real.bottom;
                            temp->area.bottom = curr->size.y - real.top;
                        }

                        temp->area.left = max(temp->area.left, 0);
                        temp->area.top = max(temp->area.top, 0);
                        if (!temp->area.right) temp->area.right = curr->size.x;
                        temp->area.right = min(temp->area.right, curr->size.x);
                        if (!temp->area.bottom) temp->area.bottom = curr->size.y;
                        temp->area.bottom = min(temp->area.bottom, curr->size.y);

                        if (temp->func) {
                            temp->func(temp, FALSE);
                            CreateThread(NULL, 0, FilterThread, temp, 0, NULL);
                            return FALSE;
                        }
                        else {
                            curr = curr->prev;
                            FreePictTree(&curr->next);
                            flgPaint = TRUE;
                            InvalidateRect(mainSurface, NULL, FALSE);
                            MessageBoxIcon(hDlg, "Функция-обработчик не задана!", NULL, MB_OK | MB_USERICON, (LPSTR)ICN_MAIN);
                        }
                    }
                    else {
                        NotifyBusy(hDlg);
                        return FALSE;
                    }
                }

                case TBT_BACK:;
                    InvalidateRect(mainSurface, NULL, FALSE);
            }

        case WM_CLOSE:
            if (!flgPaint) {
                NotifyBusy(hDlg);
                return TRUE;
            }
            FreeConf(hDlg);
            EndDialog(hDlg, 0);
            return FALSE;

        default:
            return FALSE;
    }
}



BOOL CALLBACK DialogProc(HWND hDlg, UINT uMsg, WPARAM wPrm, LPARAM lPrm) {
    RECT rcc, rcw;
    HWND hwc;

    switch (uMsg) {
        case WM_CLOSE:
            flgPaint = TRUE;
            while(conf)
                SendMessage(conf->hwid, WM_CLOSE, 0, 0);
            FreePictTree(&root);
            DeleteObject(bkBrush);
            EndDialog(hDlg, 0);
            return FALSE;


        case WM_INITDIALOG:
            mainSurface = GetDlgItem(hDlg, MIB_MAIN);
            fwdButton   = GetDlgItem(hDlg, MBT_MFWD);
            bwdButton   = GetDlgItem(hDlg, MBT_MBWD);
            progBar     = GetDlgItem(hDlg, MPB_PRCS);
            root = NULL;
            idRect = NULL;
            flgPaint = FALSE;
            curStd = GetCursor();
            curSize = LoadCursor(0, IDC_CROSS);
            curMove = LoadCursor(0, IDC_SIZEALL);
            bkBrush = CreateSolidBrush(((CLR_BACK & 0xFF) << 16) | ((CLR_BACK & 0xFF0000) >> 16) | ((CLR_BACK & 0xFF00)));
            SendMessage(hDlg, WM_SETICON, ICON_BIG, (LPARAM)LoadIcon(base, MAKEINTRESOURCE(ICN_MAIN)));
            hwc = GetDlgItem(hDlg, MCB_FLTR);
            for (rcc.left = 0; rcc.left < sizeof(allFilters)/sizeof(FLTR); rcc.left++)
                SendMessage(hwc, CB_ADDSTRING, 0, (LPARAM)allFilters[rcc.left].name);

            SendMessage(hwc, CB_INSERTSTRING, 0, (LPARAM)greeting);
            SendMessage(hwc, CB_SETCURSEL, 0, 0);

            GetWindowRect(hDlg, &rcw);
            dlgMinSz.x = rcw.right - rcw.left;
            dlgMinSz.y = rcw.bottom - rcw.top;

            GetWindowRect(GetDlgItem(hDlg, MPB_PRCS), &rcc);
            pbOffset.x = rcw.right - rcc.right;
            pbOffset.y = rcc.top - rcw.top;

            GetWindowRect(mainSurface, &rcc);
            ibOffset.x = rcw.right - rcc.right;
            ibOffset.y = rcw.bottom - rcc.bottom;

            return TRUE;


        case WM_SIZING:
            #define rcm ((RECT*)lPrm)
            if (rcm->right - rcm->left < dlgMinSz.x) {
                if (wPrm == WMSZ_LEFT || wPrm == WMSZ_TOPLEFT || wPrm == WMSZ_BOTTOMLEFT)
                    rcm->left = rcm->right - dlgMinSz.x;
                else
                    rcm->right = rcm->left + dlgMinSz.x;
            }
            if (rcm->bottom - rcm->top < dlgMinSz.y) {
                if (wPrm == WMSZ_TOP || wPrm == WMSZ_TOPLEFT || wPrm == WMSZ_TOPRIGHT)
                    rcm->top = rcm->bottom - dlgMinSz.y;
                else
                    rcm->bottom = rcm->top + dlgMinSz.y;
            }
            #undef rcm


        case WM_SIZE:
            GetWindowRect(hDlg, &rcw);

            GetWindowRect(hwc = GetDlgItem(hDlg, MPB_PRCS), &rcc);
            rcc.right = rcw.right - pbOffset.x - rcc.left;
            rcc.bottom -= rcc.top;
            ScreenToClient(hDlg, (LPPOINT)&rcc);
            MoveWindow(hwc, rcc.left, rcc.top, rcc.right, rcc.bottom, TRUE);

            GetWindowRect(mainSurface, &rcc);
            rcc.right = rcw.right - ibOffset.x - rcc.left;
            rcc.bottom = rcw.bottom - ibOffset.y - rcc.top;
            ScreenToClient(hDlg, (LPPOINT)&rcc);
            MoveWindow(mainSurface, rcc.left, rcc.top, rcc.right, rcc.bottom, TRUE);

            return (uMsg == WM_SIZING);


        case WM_DRAWITEM:
            #define dis ((DRAWITEMSTRUCT*)lPrm)
            if (flgPaint) {
                if (cpos.x < dis->rcItem.left) {
                    rcc = dis->rcItem;
                    rcc.right = rcc.left - cpos.x;
                    FillRect(dis->hDC, &rcc, bkBrush);
                }
                if (cpos.y < dis->rcItem.top) {
                    rcc = dis->rcItem;
                    rcc.bottom = rcc.top - cpos.y;
                    FillRect(dis->hDC, &rcc, bkBrush);
                }
                if (curr->size.x - cpos.x < dis->rcItem.right) {
                    rcc = dis->rcItem;
                    rcc.left = curr->size.x - cpos.x;
                    FillRect(dis->hDC, &rcc, bkBrush);
                }
                if (curr->size.y - cpos.y < dis->rcItem.bottom) {
                    rcc = dis->rcItem;
                    rcc.top = curr->size.y - cpos.y;
                    FillRect(dis->hDC, &rcc, bkBrush);
                }
                BitBlt(dis->hDC, dis->rcItem.left, dis->rcItem.top, dis->rcItem.right, dis->rcItem.bottom, curr->devc, cpos.x, cpos.y, SRCCOPY);
            }
            else FillRect(dis->hDC, &dis->rcItem, bkBrush);

            CONF *temp = conf;
            while (temp) {
                if (temp->area.left != temp->area.right && temp->area.top != temp->area.bottom) {
                    SelectObject(dis->hDC, temp->mark);

                    Rectangle(dis->hDC, temp->area.left  - cpos.x,            temp->area.top    - cpos.y - LIN_SIZE, temp->area.right - cpos.x,            temp->area.top    - cpos.y + LIN_SIZE);
                    Rectangle(dis->hDC, temp->area.left  - cpos.x,            temp->area.bottom - cpos.y - LIN_SIZE, temp->area.right - cpos.x,            temp->area.bottom - cpos.y + LIN_SIZE);
                    Rectangle(dis->hDC, temp->area.left  - cpos.x - LIN_SIZE, temp->area.top    - cpos.y,            temp->area.left  - cpos.x + LIN_SIZE, temp->area.bottom - cpos.y);
                    Rectangle(dis->hDC, temp->area.right - cpos.x - LIN_SIZE, temp->area.top    - cpos.y,            temp->area.right - cpos.x + LIN_SIZE, temp->area.bottom - cpos.y);

                    Rectangle(dis->hDC, temp->area.left  - cpos.x - BTN_SIZE, temp->area.top    - cpos.y - BTN_SIZE, temp->area.left  - cpos.x + BTN_SIZE, temp->area.top    - cpos.y + BTN_SIZE);
                    Rectangle(dis->hDC, temp->area.left  - cpos.x - BTN_SIZE, temp->area.bottom - cpos.y - BTN_SIZE, temp->area.left  - cpos.x + BTN_SIZE, temp->area.bottom - cpos.y + BTN_SIZE);
                    Rectangle(dis->hDC, temp->area.right - cpos.x - BTN_SIZE, temp->area.top    - cpos.y - BTN_SIZE, temp->area.right - cpos.x + BTN_SIZE, temp->area.top    - cpos.y + BTN_SIZE);
                    Rectangle(dis->hDC, temp->area.right - cpos.x - BTN_SIZE, temp->area.bottom - cpos.y - BTN_SIZE, temp->area.right - cpos.x + BTN_SIZE, temp->area.bottom - cpos.y + BTN_SIZE);
                }
                temp = temp->next;
            }

            return TRUE;
            #undef dis


        case WM_LBUTTONUP:
            if (idRect) {
                ReleaseCapture();
                SendMessage(GetDlgItem(idRect, TCB_RECT), BM_SETCHECK, BST_UNCHECKED, 0);
                idRect = NULL;
                SetCursor(curStd);
                InvalidateRect(mainSurface, NULL, FALSE);
            }
            return FALSE;


        case WM_MOUSEMOVE:
            if (idRect && (GetCapture() == hDlg)) {
                POINT fpos;
                GetCursorPos(&fpos);
                if (idRect == (HWND)-1) {
                    cpos.x = ppos.x - fpos.x;
                    cpos.y = ppos.y - fpos.y;
                }
                else {
                    CONF *temp = FindConf(idRect);
                    ScreenToClient(mainSurface, &fpos);
                    temp->area.right  = fpos.x + cpos.x;
                    temp->area.bottom = fpos.y + cpos.y;
                }
                InvalidateRect(mainSurface, NULL, FALSE);
            }
            return FALSE;


        case WM_COMMAND:
            switch (LOWORD(wPrm)) {
                case MIB_MAIN:
                    if (root && flgPaint) {
                        GetCursorPos(&ppos);
                        ppos.x += cpos.x;
                        ppos.y += cpos.y;
                        if (!idRect) {
                            idRect = (HWND)-1;
                            SetCursor(curMove);
                        }
                        else if (idRect != (HWND)-1) {
                            CONF *temp = FindConf(idRect);
                            ScreenToClient(mainSurface, &ppos);
                            temp->area.right  = temp->area.left = ppos.x;
                            temp->area.bottom = temp->area.top  = ppos.y;
                            SetCursor(curSize);
                        }
                        SetCapture(hDlg);
                    }
                    return FALSE;

                case MBT_DONE:
                    return DialogBoxParam(base, MAKEINTRESOURCE(DLG_CONF), NULL, (DLGPROC)ConfProc, SendMessage(GetDlgItem(hDlg, MCB_FLTR), CB_GETCURSEL, 0, 0));

                case MBT_MBWD:
                    if (root && !flgPaint) {
                        NotifyBusy(hDlg);
                        return FALSE;
                    }
                    EnableWindow(fwdButton, TRUE);
                    if (curr->prev)
                        curr = curr->prev;
                    if (!curr->prev) {
                        EnableWindow(bwdButton, FALSE);
                        SetFocus(fwdButton);
                    }
                    InvalidateRect(mainSurface, NULL, FALSE);
                    return FALSE;

                case MBT_MFWD:
                    if (root && !flgPaint) {
                        NotifyBusy(hDlg);
                        return FALSE;
                    }
                    EnableWindow(bwdButton, TRUE);
                    if (curr->next)
                        curr = curr->next;
                    if (!curr->next) {
                        EnableWindow(fwdButton, FALSE);
                        SetFocus(bwdButton);
                    }
                    InvalidateRect(mainSurface, NULL, FALSE);
                    return FALSE;

                case MBT_OPEN: {
                    if (root && !flgPaint) {
                        NotifyBusy(hDlg);
                        return FALSE;
                    }
                    char fnm[MAX_PATH + 1] = {};
                    OPENFILENAME ofn = {};

                    ofn.lStructSize = sizeof(OPENFILENAME);
                    ofn.hwndOwner = hDlg;
                    ofn.hInstance = base;
                    ofn.lpstrFilter = "BMP\0*.bmp\0GIF\0*.gif\0JPG\0*.jpg\0Все файлы\0*.*\0\0";
                    ofn.nFilterIndex = 4;
                    ofn.lpstrFile = (LPSTR)&fnm;
                    ofn.nMaxFile = MAX_PATH;
                    ofn.lpstrTitle = "Загрузка картинки";
                    ofn.Flags = OFN_ENABLESIZING | OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST | OFN_HIDEREADONLY;

                    if (GetOpenFileName(&ofn))
                        LoadRootPict(fnm, hDlg);
                    return FALSE;
                }

                case MBT_SAVE: {
                    if (root && !flgPaint) {
                        NotifyBusy(hDlg);
                        return FALSE;
                    }
                    char fnm[MAX_PATH + 1] = {};
                    OPENFILENAME ofn = {};

                    ofn.lStructSize = sizeof(OPENFILENAME);
                    ofn.hwndOwner = hDlg;
                    ofn.hInstance = base;
                    ofn.lpstrFilter = "BMP\0*.bmp\0Все файлы\0*.*\0\0";
                    ofn.nFilterIndex = 0;
                    ofn.lpstrFile = (LPSTR)&fnm;
                    ofn.nMaxFile = MAX_PATH;
                    ofn.lpstrTitle = "Сохранение картинки";
                    ofn.lpstrDefExt = "bmp";
                    ofn.Flags = OFN_EXTENSIONDIFFERENT | OFN_ENABLESIZING | OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST | OFN_HIDEREADONLY;

                    if (GetSaveFileName(&ofn))
                        SaveCurrPict(fnm, hDlg);
                    return FALSE;
                }
            }


        default:
            return FALSE;
    }
}



int CALLBACK WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nShowCmd) {
    SYSTEM_INFO syin;
    GetSystemInfo(&syin);
    numFreeCores = syin.dwNumberOfProcessors;
    if (!numFreeCores) numFreeCores = 1;
    if (!(nShowCmd = strlen(lpCmdLine))) {
        INITCOMMONCONTROLSEX icc;
        base = hInstance;

        srand(time(0));
        icc.dwSize = sizeof(INITCOMMONCONTROLSEX);
        icc.dwICC = ICC_STANDARD_CLASSES | ICC_LISTVIEW_CLASSES | ICC_PROGRESS_CLASS | ICC_UPDOWN_CLASS;
        InitCommonControlsEx(&icc);

        return DialogBoxParam(base, MAKEINTRESOURCE(DLG_MAIN), NULL, (DLGPROC)DialogProc, 0);
    }

    LPSTR pcmd = strdup(lpCmdLine);
    LONG fout, finp, fctl;
    CONF vcnf = {};
    FLOAT coef;

    for (fout = nShowCmd;   (fout > 0) && (pcmd[fout - 1] == ' '); fout--);
    for (pcmd[fout] = '\0'; (fout > 0) && (pcmd[fout - 1] != ' '); fout--);
    for (finp = fout;       (finp > 0) && (pcmd[finp - 1] == ' '); finp--);
    for (pcmd[finp] = '\0'; (finp > 0) && (pcmd[finp - 1] != ' '); finp--);
    if (finp < 3) return -1;

    for (fctl = 0; (fctl < finp) && (pcmd[fctl] != '-'); fctl++);
    if (finp - fctl < 3) return -2;
    pcmd[finp - 1] = '\0';

    if (!LoadRootPict(pcmd + finp, NULL)) return -3;
    finp = pcmd[fctl + 1];
    fctl += 2;
    while (pcmd[fctl] == ' ') fctl++;

    if (finp - fctl > 1)  {
        DuplicatePict(curr);
        curr = curr->next;
        vcnf.area.right = curr->size.x;
        vcnf.area.bottom = curr->size.y;

        switch (finp) {
            case 'g':
                sscanf(pcmd + fctl, "%f", &coef);
                vcnf.data = tr(fabs(coef * 1000.0));
                vcnf.data %= 100001;
                GaussBlur(&vcnf, TRUE);
                break;

            case 'm':
                sscanf(pcmd + fctl, "%f", &coef);
                vcnf.data = tr(fabs(coef));
                vcnf.data &= 0xFF;
                MedianFlt(&vcnf, TRUE);
                break;

            case 's':
                vcnf.data = 2;
                SobelEdge(&vcnf, TRUE);
                break;

            case 'a':
                GrayWorld(&vcnf, TRUE);
                break;

            case 'x':
                vcnf.data = 1;
                FadToGray(&vcnf, TRUE);
                break;

            case 'c':
                vcnf.data = 5 | 256;
                AContLvls(&vcnf, TRUE);
                break;

            case 'l':
                vcnf.data = 5 | 512;
                AContLvls(&vcnf, TRUE);
                break;

            case 'r':
                sscanf(pcmd + fctl, "%f", &coef);
                vcnf.data = tr(fabs(coef));
                vcnf.data &= 0xFF;
                vcnf.data |= 0x200;
                if (coef < 0) vcnf.data |= 0x100;
                RotateFlt(&vcnf, TRUE);
                break;

            case 'z':
                sscanf(pcmd + fctl, "%f", &coef);
                vcnf.data = tr(fabs(coef * 1000.0));
                vcnf.data &= 0xFFFF;
                vcnf.data |= 0x20000;
                ScalerFlt(&vcnf, TRUE);
                break;

            case 'k':
                for (finp = fctl + 1; (finp < fout - 1) && (pcmd[finp] != '\''); finp++);
                pcmd[finp] = '\0';
                vcnf.data = (DWORD)pcmd + fctl + 1;
                CustomLin(&vcnf, TRUE);
                break;

            case 'f':
                sscanf(pcmd + fctl, "%f", &coef);
                vcnf.data = tr(fabs(coef));
                vcnf.data %= 101;
                BrkGlsFlt(&vcnf, TRUE);
                break;

            case 'h':
                sscanf(pcmd + fctl, "%f", &coef);
                vcnf.data = tr(fabs(coef));
                vcnf.data %= 201;
                vcnf.data |= 0x100;
                SineWaves(&vcnf, TRUE);
                break;

            case 'w':
                sscanf(pcmd + fctl, "%f", &coef);
                vcnf.data = tr(fabs(coef));
                vcnf.data %= 201;
                vcnf.data |= 0x200;
                SineWaves(&vcnf, TRUE);
                break;

            case 't':
                sscanf(pcmd + fctl, "%f", &coef);
                vcnf.data = tr(fabs(coef));
                vcnf.data %= 101;
                GlassTile(&vcnf, TRUE);
                break;

            default:
                finp = -4;
        }
        if (finp >= 0) SaveCurrPict(pcmd + fout, NULL);
    }
    else finp = -2;

    FreePictTree(&root);
    return finp;
}
