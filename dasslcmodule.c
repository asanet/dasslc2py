/*##################################################################
#                                                                  #
#   This is the source code of the python wrapper for the          #
#   "Differential-Algebraic System Solver  in C" by A. R. Secchi   #
#                                                                  #
#   Author: Ataide Neto                                            #
#   email: ataide@peq.coppe.ufrj.br                                #
#   Universidade Federal do Rio de Janeiro                         #
#   Version: 0.1-5                                                 #
#                                                                  #
##################################################################*/

/* Change log: 
 *
 * v0.1-0 First working version in python2
 * v0.1-1 Added: python3 compatibility
 * v0.1-2 Added: inputfile and user-defined jacobian support
 * v0.1-3 Added: steady state and sparse algebra support
 * v0.1-4 Added: user-defined jacobian support (sparse)
 * v0.1-5 Fixed: sparse algebra compilation
 */

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <numpy/arrayobject.h>
#include "dasslc/dasslc.h"

#define FREEALL() Py_XDECREF(t_array); Py_XDECREF(y_array);\
                  Py_XDECREF(yp_array); Py_XDECREF(idx_array);\
                  if (y) free(y);\
                  if (yp) free(yp);\
                  if (index) free(index);\
                  daFree(&root);

#define REQS NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED|NPY_ARRAY_WRITEABLE|\
             NPY_ARRAY_ENSUREARRAY|NPY_ARRAY_ENSURECOPY


// The function's prototypes */
static PyObject* dasslc_solve(PyObject *self, PyObject *args);

// The method's table
static PyMethodDef dasslcMethods[] = {
    {"solve", dasslc_solve, METH_VARARGS , "Solve the problem."},
    {NULL, NULL, 0, NULL}
};

// The  global variables
static DASSLC_RES residuals;         //The C residual function
static DASSLC_JAC jacobian;          //The C jacobian function
static PyObject *pyres = NULL;       //The Python residual function
static PyObject *pyjac = NULL;       //The Python jacobian function


#if PY_MAJOR_VERSION >=3
    // The module definition structure PY3
    static struct PyModuleDef dasslcmodule = {
        PyModuleDef_HEAD_INIT,
        "dasslc",     // module name
        NULL,         // module documentation
        -1,           // size of per-interpreter state of the module,
                      // or -1 if the module keeps state in global variables
        dasslcMethods //Methods
    };
    // The module initialization function PY3
    PyMODINIT_FUNC PyInit_dasslc(void){

        PyObject *m = PyModule_Create(&dasslcmodule);
        if (m == NULL)
            return NULL;

        import_array();
        return m;
    }
#else 
    // The module initialization function PY2
    PyMODINIT_FUNC initdasslc(void){
        PyObject *m = Py_InitModule("dasslc", dasslcMethods);
        if (m == NULL)
            return;

        // Load numpy functionality.
        import_array();
    }
#endif

// The function declaration
static PyObject* dasslc_solve(PyObject *self, PyObject *args){
    //python call: dasslc.solve(resfun, tspan, y0, yp0, rpar, rtol, atol, index, inputfile, jac)

    // The memory allocation
    int ntp = -1, ntp2 = -1, neq = -1, ndr = -1, idxnum = -1;
    int *index = NULL, i = 0, j = 0;
    double t0 = 0, tf = 0, dt = 0, *y = NULL, *yp = NULL;
    PyArrayObject *t_sol = NULL, *y_sol = NULL, *yp_sol = NULL;
    PyArrayObject *t_array = NULL, *y_array = NULL, *yp_array = NULL;
    PyArrayObject *idx_array = NULL;
    PyObject *result = NULL, *arglist = NULL;
    char *inputfile = NULL;
    PTR_ROOT root;
    BOOL err = 0;

    // The python inputs and outputs
    // >> mandatory:
    PyObject *t_obj = NULL, *y_obj = NULL, *resfun_obj = NULL; 
    // >> optional:
    PyObject *yp_obj = NULL, *rpar_obj = NULL, *idx_obj = NULL, *jac_obj = NULL;
    double atol = 1e-10, rtol = 1e-8; 

    // Parse inputs
    if (!PyArg_ParseTuple(args, "OOO|OOddOzO", &resfun_obj, &t_obj, &y_obj,
                          &yp_obj, &rpar_obj, &rtol, &atol, &idx_obj, &inputfile, &jac_obj))
        return NULL;

    // Interpret the input objects as numpy arrays.
    t_array = (PyArrayObject*)PyArray_FROMANY(t_obj, NPY_DOUBLE, 0, 1, REQS);
    y_array = (PyArrayObject*)PyArray_FROMANY(y_obj, NPY_DOUBLE, 0, 1, REQS);

    if (yp_obj && yp_obj != Py_None)
        yp_array = (PyArrayObject*)PyArray_FROMANY(yp_obj, NPY_DOUBLE, 0, 1, REQS);
    if (idx_obj && idx_obj != Py_None)
        idx_array = (PyArrayObject*)PyArray_FROMANY(idx_obj, NPY_DOUBLE, 0, 1, REQS);

    // If that didn't work, throw an exception.
    if (t_array == NULL || y_array == NULL){
        FREEALL();
        PyErr_SetString(PyExc_TypeError, "t and y must be 1D-arrays");
        return NULL;
    }

    // Get dimensions
    ntp = PyArray_NDIM(t_array) == 0 ? 1 : (int)PyArray_DIM(t_array, 0);
    neq = PyArray_NDIM(y_array) == 0 ? 1 : (int)PyArray_DIM(y_array, 0);

    // Check if residual function is callable
    if (!PyCallable_Check(resfun_obj)){
        FREEALL();
        PyErr_SetString(PyExc_TypeError, "Cannot call provided residual function.");
        return NULL;
    }
    Py_XINCREF(resfun_obj);      //Add a reference to new callback
    Py_XDECREF(pyres);           //Dispose of previous callback
    pyres = resfun_obj;          //Remember new callback

    if (rpar_obj && rpar_obj != Py_None)
        arglist = Py_BuildValue("dOOO",3.14,y_obj,y_obj,rpar_obj);
    else
        arglist = Py_BuildValue("dOO",3.14,y_obj,y_obj);

    result = PyObject_CallObject(pyres, arglist);
    PyObject *dummyO;
    int dummyi;
    if (!result || !PyArg_ParseTuple(result, "Oi",&dummyO,&dummyi)){
        FREEALL();
        Py_XDECREF(result);
        Py_XDECREF(arglist);
        Py_XDECREF(dummyO);
        PyErr_SetString(PyExc_TypeError, "There's some problem in residual function. Check I/O!");
        return NULL;
    }
    PyArrayObject *dummyVec = (PyArrayObject*)PyArray_FROMANY(dummyO, NPY_DOUBLE, 0, 1, REQS);
    int dummyVar = PyArray_NDIM(dummyVec) == 0 ? 1 : (int)PyArray_DIM(dummyVec, 0);
    if (dummyVar != neq){
        FREEALL();
        Py_XDECREF(result);
        Py_XDECREF(arglist);
        Py_XDECREF(dummyVec);
        Py_XDECREF(dummyO);
        PyErr_SetString(PyExc_TypeError, "Residual function must return a vector with the same length of y0!");
        return NULL;
    }

    Py_XDECREF(dummyO);
    Py_XDECREF(dummyVec);
    Py_XDECREF(arglist);
    //Py_XDECREF(result);

    // Check if jacobian function is callable
    if (jac_obj && jac_obj != Py_None){
        if (!PyCallable_Check(jac_obj)){
            FREEALL();
            PyErr_SetString(PyExc_TypeError, "Cannot call provided jacobian function.");
            return NULL;
        }
        Py_XINCREF(jac_obj);
        Py_XDECREF(pyjac);
        pyjac = jac_obj;

        if (rpar_obj && rpar_obj != Py_None)
            arglist = Py_BuildValue("dOOdO",3.14,y_obj,y_obj,2.71,rpar_obj);
        else
            arglist = Py_BuildValue("dOOd",3.14,y_obj,y_obj,2.71);

        result = PyObject_CallObject(pyjac, arglist);
        PyObject *dummyI = NULL, *dummyJ = NULL;
        if (!result || !PyArg_ParseTuple(result, "Oi|OO",&dummyO,&dummyi,&dummyI,&dummyJ)){
            FREEALL();
            Py_XDECREF(result);
            Py_XDECREF(arglist);
            PyErr_SetString(PyExc_TypeError, "There's some problem in jacobian function. Check I/O!");
            return NULL;
        }

        PyArrayObject *dummyMat = (PyArrayObject*)PyArray_FROMANY(dummyO, NPY_DOUBLE, 0, 2, REQS);
        int dummyN = PyArray_NDIM(dummyMat) == 0 ? 1 : (int)PyArray_DIM(dummyMat, 0);
        int dummyM = PyArray_NDIM(dummyMat) < 2  ? 1 : (int)PyArray_DIM(dummyMat, 1);
        if (dummyN*dummyM != neq*neq || PyArray_NDIM(dummyMat) == 1){
            FREEALL();
            Py_XDECREF(result);
            Py_XDECREF(arglist);
            Py_XDECREF(dummyMat);
            Py_XDECREF(dummyO);
            Py_XDECREF(dummyI);
            Py_XDECREF(dummyJ);
            PyErr_SetString(PyExc_TypeError, "Jacobian function must return a neq-by-neq matrix!");
            return NULL;
        }
        //Py_XDECREF(result);
        Py_XDECREF(arglist);
        Py_XDECREF(dummyMat);
        Py_XDECREF(dummyO);
        Py_XDECREF(dummyI);
        Py_XDECREF(dummyJ);
    }

    // Get pointers to the data as C-types.
    y = (double*) malloc(neq*sizeof(double));
    for (j = 0; j < neq; j++)
        y[j] = *(double*)PyArray_GETPTR1(y_array,j);

    if (yp_array){
        ndr = PyArray_NDIM(yp_array) == 0 ? 1 : (int)PyArray_DIM(yp_array, 0);
        if (ndr >= 0 && ndr != neq){ //Throw an exception
            FREEALL();
            PyErr_SetString(PyExc_TypeError, "yp0 must have the same length of y0!");
            return NULL;
        }
        yp = (double*) malloc(neq*sizeof(double));
        for (j = 0; j < neq; j++)
            yp[j] = *(double*)PyArray_GETPTR1(yp_array,j);
    }
    if (idx_array){
        idxnum = PyArray_NDIM(idx_array) == 0 ? 1 : (int)PyArray_DIM(idx_array, 0);
        if (idxnum >= 0 && idxnum != neq){//Throw an exception
            FREEALL();
            PyErr_SetString(PyExc_TypeError, "Index vector must have the same length of y0!");
            return NULL;
        }
        index = (int*) malloc(neq*sizeof(int));
        for (j = 0; j < neq; j++)
            index[j] = (int) *(double*)PyArray_GETPTR1(idx_array,j);
            //Workaround: Problem creating an int array directly in python
            //            So, create a double array then convert it to int here
    }

    // Set the rpar if any
    if (!rpar_obj || rpar_obj == Py_None)
        root.user = NULL;
    else
        root.user = (void*)rpar_obj;

    // Create the solution vector
    ntp2 = ntp > 2 ? ntp : 100;
    npy_intp dims1[1] = {ntp2};
    npy_intp dims2[2] = {ntp2,neq};
    t_sol = (PyArrayObject*)PyArray_EMPTY(1,dims1,NPY_DOUBLE,0);
    y_sol = (PyArrayObject*)PyArray_EMPTY(2,dims2,NPY_DOUBLE,0);
    yp_sol = (PyArrayObject*)PyArray_EMPTY(2,dims2,NPY_DOUBLE,0);

    // Check input file
    if (inputfile && !fopen(inputfile,"r")){
        FREEALL();
        Py_XDECREF(t_sol);
        Py_XDECREF(y_sol);
        Py_XDECREF(yp_sol);
        PyErr_SetString(PyExc_TypeError, "Cannot open provided inputfile!");
        return NULL;
    }else if (!inputfile){
        inputfile = "?";
        if (jac_obj)
            PyErr_WarnEx(PyExc_RuntimeWarning,
                         "Ignoring provided jacobian because inputfile is not properly configured.",
                         1);
    }

    // Call the daSetup function
    t0 = ntp == 1 ? 0 : *(double*)PyArray_GETPTR1(t_array,0);
    err = daSetup(inputfile,&root,residuals,neq,t0,y,yp,index,NULL,NULL);
    if (err){
        FREEALL();
        Py_XDECREF(t_sol);
        Py_XDECREF(y_sol);
        Py_XDECREF(yp_sol);
        char buff[128] = "Setup error: ";
        sprintf(buff,"%d",err);
        PyErr_SetString(PyExc_TypeError, buff);
        return NULL;
    }

    // Configure root structure
    root.iter.stol = 1;
    root.iter.atol[0] = atol;
    root.iter.rtol[0] = rtol;

    // Find initial derivatives if not given
    if (ntp == 1 && t_obj != Py_None){
        dt = (double) *(double*)PyArray_GETPTR1(t_array,0)/(ntp2-1);
        tf = t0 + dt;
    }else if (ntp == 2){
        dt = (double) (*(double*)PyArray_GETPTR1(t_array,1) - *(double*)PyArray_GETPTR1(t_array,0))/(ntp2-1);
        tf = t0 + dt;
    }else if (t_obj == Py_None){
        tf = 1;
    }else{
        tf = *(double*)PyArray_GETPTR1(t_array,1);
    }

    if (yp == NULL && t_obj != Py_None ){
        err = dasslc(INITIAL_COND, &root, residuals, &t0, tf, pyjac ? jacobian : NULL, NULL);
        if (err < 0){
            FREEALL();
            Py_XDECREF(t_sol);
            Py_XDECREF(y_sol);
            Py_XDECREF(yp_sol);
            char buff[128] = "Failed in finding consistent initial condition. Error: ";
            sprintf(buff,"%d",err);
            PyErr_SetString(PyExc_TypeError, buff);
            return NULL;
        }
    }

    // Update soluton vector
    if (t_obj != Py_None){
        *(double*)PyArray_GETPTR1(t_sol,0) = root.t;
        for (j = 0; j < neq; j++){
            *(double*)PyArray_GETPTR2(y_sol,0,j) = root.y[j];
            *(double*)PyArray_GETPTR2(yp_sol,0,j) = root.yp[j];
        }

        // Call the dasslc function for all tspan
        for (i = 1; i < ntp2; i++){
            tf = ntp > 2 ? *(double*)PyArray_GETPTR1(t_array,i) : t0 + dt;
            err = dasslc(TRANSIENT, &root, residuals, &t0, tf, pyjac ? jacobian : NULL, NULL);
            if (err < 0){
                FREEALL();
                Py_XDECREF(t_sol);
                Py_XDECREF(y_sol);
                Py_XDECREF(yp_sol);
                char buff[128] = "Error during integration: ";
                sprintf(buff,"%d",err);
                PyErr_SetString(PyExc_TypeError, buff);
                return NULL;
            }
            *(double*)PyArray_GETPTR1(t_sol,i) = root.t;
            for (j = 0; j < neq; j++){
                *(double*)PyArray_GETPTR2(y_sol,i,j) = root.y[j];
                *(double*)PyArray_GETPTR2(yp_sol,i,j) = root.yp[j];
            }
        }
    }else{
        dims1[0] = 1;
        npy_intp dims0[1] = {neq};
        t_sol = (PyArrayObject*)PyArray_EMPTY(1,dims1,NPY_DOUBLE,0);
        y_sol = (PyArrayObject*)PyArray_EMPTY(1,dims0,NPY_DOUBLE,0);
        yp_sol = (PyArrayObject*)PyArray_EMPTY(1,dims0,NPY_DOUBLE,0);
        err = dasslc(STEADY_STATE, &root, residuals, &t0, tf, pyjac ? jacobian : NULL, NULL);
        if (err < 0){
            FREEALL();
            Py_XDECREF(t_sol);
            Py_XDECREF(y_sol);
            Py_XDECREF(yp_sol);
            char buff[128] = "Error in finding steady state: ";
            sprintf(buff,"%d",err);
            PyErr_SetString(PyExc_TypeError, buff);
            return NULL;
        }
        *(double*)PyArray_GETPTR1(t_sol,0) = 0;
        for (j = 0; j < neq; j++){
            *(double*)PyArray_GETPTR1(y_sol,j) = root.y[j];
            *(double*)PyArray_GETPTR1(yp_sol,j) = root.yp[j];
        }
    }

    // Print dasslc output log
    daStat (root.savefile, &root);

    // Clean Up
    FREEALL();

    // Build the output tuple
    return Py_BuildValue("NNN", t_sol, y_sol, yp_sol);
}

static BOOL residuals(PTR_ROOT *root, REAL t, REAL *y, REAL *yp, REAL *res, BOOL *jac){
    //Interface with the python residual function

    // Memory allocation
    PyObject *arglist = NULL, *result = NULL, *res_obj = NULL;
    PyArrayObject *res_array = NULL, *y_array = NULL, *yp_array = NULL;
    int ires = -1, i = 0;

    // Build the arglist (convert c-array to PyArray)
    int neq = root -> rank;
    npy_intp dims[1] = {neq};
    y_array = (PyArrayObject*)PyArray_EMPTY(1,dims,NPY_DOUBLE,0);
    yp_array = (PyArrayObject*)PyArray_EMPTY(1,dims,NPY_DOUBLE,0);

    for (i = 0; i < neq; i++){
        *(double*)PyArray_GETPTR1(y_array,i) = y[i];
        *(double*)PyArray_GETPTR1(yp_array,i) = yp[i];
    }

    // Parse arglist checking if rpar exists
    if (root->user)
        arglist = Py_BuildValue("dOOO",t,y_array,yp_array,(PyObject*)root->user);
    else
        arglist = Py_BuildValue("dOO",t,y_array,yp_array);

    // Call the python function
    result = PyObject_CallObject(pyres, arglist);
    Py_XDECREF(arglist);

    // Parse the result tuple
    PyArg_ParseTuple(result, "Oi", &res_obj, &ires);
    res_array = (PyArrayObject*)PyArray_FROMANY(res_obj, NPY_DOUBLE, 0, 1, REQS);

    if (ires){
        Py_XDECREF(y_array);
        Py_XDECREF(yp_array);
        Py_XDECREF(result);
        Py_XDECREF(res_array);
        return ires;
    }

    // Convert result to c-array res
    for (i = 0; i < neq; i++)
        res[i] = *(double*)PyArray_GETPTR1(res_array,i);

    // Clean up
    Py_XDECREF(y_array);
    Py_XDECREF(yp_array);
    Py_XDECREF(result);
    Py_XDECREF(res_array);

    return ires;
}

#define PD(i,j) (*(pd + neq * (i) + j))

static BOOL jacobian(PTR_ROOT *root, REAL t, REAL *y, REAL *yp, REAL cj, void *ja, DASSLC_RES *residuals){

    // Memory allocation
    PyObject *arglist = NULL, *result = NULL, *pd_obj = NULL, *i_obj = NULL, *j_obj = NULL;
    PyArrayObject *pd_array = NULL, *y_array = NULL, *yp_array = NULL, *i_list = NULL, *j_list = NULL;
    int ires = -1, i = 0, j = 0, k = 0;

    // Build the arglist (convert c-array to PyArray)
    int neq = root -> rank;
    npy_intp dims[1] = {neq};
    y_array = (PyArrayObject*)PyArray_EMPTY(1,dims,NPY_DOUBLE,0);
    yp_array = (PyArrayObject*)PyArray_EMPTY(1,dims,NPY_DOUBLE,0);

    for (i = 0; i < neq; i++){
        *(double*)PyArray_GETPTR1(y_array,i) = y[i];
        *(double*)PyArray_GETPTR1(yp_array,i) = yp[i];
    }

    // Parse arglist checking if rpar exists
    if (root->user)
        arglist = Py_BuildValue("dOOdO",t,y_array,yp_array,cj,(PyObject*)root->user);
    else
        arglist = Py_BuildValue("dOOd",t,y_array,yp_array,cj);

    // Call the python function
    result = PyObject_CallObject(pyjac, arglist);
    Py_XDECREF(arglist);

    PyArg_ParseTuple(result, "Oi|OO", &pd_obj, &ires, &i_obj, &j_obj);
    pd_array = (PyArrayObject*)PyArray_FROMANY(pd_obj, NPY_DOUBLE, 0, 2, REQS);

    if (ires){
        Py_XDECREF(y_array);
        Py_XDECREF(yp_array);
        Py_XDECREF(result);
        Py_XDECREF(pd_array);
        return ires;
    }

    // Convert result to c-array res
    switch (root->jac.mtype){
        case USER_DENSE:{
            REAL *pd = (REAL *)ja;
            for (i = 0; i < neq; i++)
                for (j = 0; j < neq; j++)
                    PD(i,j) = *(double*)PyArray_GETPTR2(pd_array,i,j);
            break;
        }
        case USER_BAND:{
            REAL *pd = (REAL *)ja;
            int m = root->jac.lband + root->jac.uband, k = 0;
            for (i = 0; i < neq; i++)
                for (j = 0; j < neq; j++){
                    k = i - j + m;
                    PD(k,j) = *(double*)PyArray_GETPTR2(pd_array,i,j);
                }
            break;
        }
#ifdef SPARSE
        case USER_SPARSE:{
            if (!i_obj || !j_obj){
                ires = 1;
                break;
            }
            i_list = (PyArrayObject*)PyArray_FROMANY(i_obj, NPY_DOUBLE, 0, 0, REQS);
            j_list = (PyArrayObject*)PyArray_FROMANY(j_obj, NPY_DOUBLE, 0, 0, REQS);
            int iN = (int)PyArray_DIM(i_list, 0);
            int jN = (int)PyArray_DIM(j_list, 0);
            if (iN != jN){
                ires = 1;
                Py_XDECREF(i_list);
                Py_XDECREF(j_list);
                break;
            }
            char *pd = (char*)ja;
            for (k = 0; k < iN; k++){
                i = (int) *(double*)PyArray_GETPTR1(i_list,k);
                j = (int) *(double*)PyArray_GETPTR1(j_list,k);
                daSparse_value(pd,i,j) = *(double*)PyArray_GETPTR2(pd_array,i,j);
            }
            Py_XDECREF(i_list);
            Py_XDECREF(j_list);
            break;
        }
#endif
    }

    // Clean up
    Py_XDECREF(y_array);
    Py_XDECREF(yp_array);
    Py_XDECREF(result);
    Py_XDECREF(pd_array);

    return ires;
}
