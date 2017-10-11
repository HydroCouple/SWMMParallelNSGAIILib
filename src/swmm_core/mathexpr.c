/******************************************************************************
**  MODULE:        MATHEXPR.C
**  PROJECT:       EPA SWMM 5.1
**  DESCRIPTION:   Evaluates symbolic mathematical expression consisting
**                 of numbers, variable names, math functions & arithmetic
**                 operators.
**  AUTHORS:       L. Rossman, US EPA - NRMRL
**                 F. Shang, University of Cincinnati
**  VERSION:       5.1.008
**  LAST UPDATE:   04/01/15
******************************************************************************/
/*
**   Operand codes:
** 	   1 = (
** 	   2 = )
** 	   3 = +
** 	   4 = - (subtraction)
** 	   5 = *
** 	   6 = /
** 	   7 = number
** 	   8 = user-defined variable
** 	   9 = - (negative)
**	  10 = cos
**	  11 = sin
**	  12 = tan
**	  13 = cot
**	  14 = abs
**	  15 = sgn
**	  16 = sqrt
**	  17 = log
**	  18 = exp
**	  19 = asin
**	  20 = acos
**	  21 = atan
**	  22 = acot
**    23 = sinh
**	  24 = cosh
**	  25 = tanh
**	  26 = coth
**	  27 = log10
**    28 = step (x<=0 ? 0 : 1)
**	  31 = ^
******************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <ctype.h>
#include <stdlib.h>
//#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mathexpr.h"
#include "globals.h"

#define MAX_STACK_SIZE  1024

//  Local declarations
//--------------------
//  Structure for binary tree representation of math expression
struct TreeNode
{
    int    opcode;                // operator code
    int    ivar;                  // variable index
    double fvalue;                // numerical value
    struct TreeNode *left;        // left sub-tree of tokenized formula
    struct TreeNode *right;       // right sub-tree of tokenized formula
};
typedef struct TreeNode ExprTree;

// Local variables
//----------------
//static int    project->Err;
//static int    project->Bc;
//static int    project->PrevLex, project->CurLex;
//static int    project->Len,  project->Pos;
//static char   *project->S;
//static char    project->Token[255];
//static int    project->Ivar;
//static double project->Fvalue;

// math function names
char *MathFunc[] =  {"COS", "SIN", "TAN", "COT", "ABS", "SGN",
                     "SQRT", "LOG", "EXP", "ASIN", "ACOS", "ATAN",
                     "ACOT", "SINH", "COSH", "TANH", "COTH", "LOG10",
                     "STEP", NULL};

// Local functions
//----------------
static int        sametext(char *, char *);
static int        isDigit(char);
static int        isLetter(char);
static void       getToken(Project *project);
static int        getMathFunc(Project *project);
static int        getVariable(Project* project);
static int        getOperand(Project* project);
static int        getLex(Project* project);
static double     getNumber(Project* project);
static ExprTree * newNode(Project *project);
static ExprTree * getSingleOp(Project* project, int *);
static ExprTree * getOp(Project* project, int *);
static ExprTree * getTree(Project* project);
static void       traverseTree(ExprTree *, MathExpr **);
static void       deleteTree(ExprTree *);

// Callback functions
static int    (*getVariableIndex) (Project*, char *); // return index of named variable

//=============================================================================

int  sametext(char *s1, char *s2)
/*
**  Purpose:
**    performs case insensitive comparison of two strings.
**
**  Input:
**    s1 = character string
**    s2 = character string.
**  
**  Returns:
**    1 if strings are the same, 0 otherwise.
*/
{
   int i;
   for (i=0; toupper(s1[i]) == toupper(s2[i]); i++)
     if (!s1[i+1] && !s2[i+1]) return(1);
   return(0);
}

//=============================================================================

int isDigit(char c)
{
    if (c >= '1' && c <= '9') return 1;
    if (c == '0') return 1;
    return 0;
}

//=============================================================================

int isLetter(char c)
{
    if (c >= 'a' && c <= 'z') return 1;
    if (c >= 'A' && c <= 'Z') return 1;
    if (c == '_') return 1;
    return 0;
}

//=============================================================================

void getToken(Project *project)
{
    char c[] = " ";
    strcpy(project->Token, "");
    while (project->Pos <= project->Len &&
        ( isLetter(project->S[project->Pos]) || isDigit(project->S[project->Pos]) ) )
    {
        c[0] = project->S[ project->Pos];
        strcat( project->Token, c);
         project->Pos++;
    }
     project->Pos--;
}

//=============================================================================

int getMathFunc(Project *project)
{
    int i = 0;
    while (MathFunc[i] != NULL)
    {
        if (sametext(MathFunc[i],  project->Token)) return i+10;
        i++;
    }
    return(0);
}

//=============================================================================

int getVariable(Project* project)
{
    if ( !getVariableIndex ) return 0;
    project->Ivar = getVariableIndex(project, project->Token);
    if (project->Ivar >= 0) return 8;
    return 0;
}

//=============================================================================

double getNumber(Project* project)
{
    char c[] = " ";
    char sNumber[255];
    int  errflag = 0;

    /* --- get whole number portion of number */
    strcpy(sNumber, "");
    while ( project->Pos < project->Len && isDigit(project->S[ project->Pos]))
    {
        c[0] = project->S[ project->Pos];
        strcat(sNumber, c);
         project->Pos++;
    }

    /* --- get fractional portion of number */
    if ( project->Pos < project->Len)
    {
        if (project->S[ project->Pos] == '.')
        {
            strcat(sNumber, ".");
             project->Pos++;
            while ( project->Pos < project->Len && isDigit(project->S[ project->Pos]))
            {
                c[0] = project->S[ project->Pos];
                strcat(sNumber, c);  
                 project->Pos++;
            }
        }

        /* --- get exponent */
        if ( project->Pos < project->Len && (project->S[ project->Pos] == 'e' || project->S[ project->Pos] == 'E'))
        {
            strcat(sNumber, "E");  
             project->Pos++;
            if ( project->Pos >= project->Len) errflag = 1;
            else
            {
                if (project->S[ project->Pos] == '-' || project->S[ project->Pos] == '+')
                {
                    c[0] = project->S[ project->Pos];
                    strcat(sNumber, c);  
                     project->Pos++;
                }
                if ( project->Pos >= project->Len || !isDigit(project->S[ project->Pos])) errflag = 1;
                else while (  project->Pos < project->Len && isDigit(project->S[ project->Pos]))
                {
                    c[0] = project->S[ project->Pos];
                    strcat(sNumber, c);  
                     project->Pos++;
                }
            }
        }
    }
     project->Pos--;
    if (errflag) return 0;
    else return atof(sNumber);
}

//=============================================================================

int getOperand(Project* project)
{
    int code;
    switch(project->S[ project->Pos])
    {
      case '(': code = 1;  break;
      case ')': code = 2;  break;
      case '+': code = 3;  break;
      case '-': code = 4;
        if ( project->Pos < project->Len-1 &&
            isDigit(project->S[ project->Pos+1]) &&
            (project->CurLex == 0 || project->CurLex == 1))
        {
             project->Pos++;
            project->Fvalue = -getNumber(project);
            code = 7;
        }
        break;
      case '*': code = 5;  break;
      case '/': code = 6;  break;
      case '^': code = 31; break;
      default:  code = 0;
    }
    return code;
}

//=============================================================================

int getLex(Project* project)
{
    int n;

    /* --- skip spaces */
    while (  project->Pos < project->Len && project->S[ project->Pos] == ' ' )  project->Pos++;
    if (  project->Pos >= project->Len ) return 0;

    /* --- check for operand */
    n = getOperand(project);

    /* --- check for function/variable/number */
    if ( n == 0 )
    {
        if ( isLetter(project->S[ project->Pos]) )
        {
            getToken(project);
            n = getMathFunc(project);
            if ( n == 0 ) n = getVariable(project);
        }
        else if ( isDigit(project->S[ project->Pos]) )
        {
            n = 7;
            project->Fvalue = getNumber(project);
        }
    }
     project->Pos++;
    project->PrevLex = project->CurLex;
    project->CurLex = n;
    return n;
}

//=============================================================================

ExprTree * newNode(Project *project)
{
    ExprTree *node;
    node = (ExprTree *) malloc(sizeof(ExprTree));
    if (!node) project->Err = 2;
    else
    {
        node->opcode = 0;
        node->ivar   = -1;
        node->fvalue = 0.;
        node->left   = NULL;
        node->right  = NULL;
    }
    return node;
}

//=============================================================================

ExprTree * getSingleOp(Project* project,int *lex)
{
    int bracket;
    int opcode;
    ExprTree *left;
    ExprTree *right;
    ExprTree *node;

    /* --- open parenthesis, so continue to grow the tree */
    if ( *lex == 1 )
    {
        project->Bc++;
		left = getTree(project);
    }

    else
    {
        /* --- Error if not a singleton operand */
        if ( *lex < 7 || *lex == 9 || *lex > 30)
        {
            project->Err = 1;
            return NULL;
        }

        opcode = *lex;

        /* --- simple number or variable name */
        if ( *lex == 7 || *lex == 8 )
        {
            left = newNode(project);
            left->opcode = opcode;
            if ( *lex == 7 ) left->fvalue = project->Fvalue;
            if ( *lex == 8 ) left->ivar = project->Ivar;
        }

        /* --- function which must have a '(' after it */
        else
        {
            *lex = getLex(project);
            if ( *lex != 1 )
            {
               project->Err = 1;
               return NULL;
            }
            project->Bc++;
            left = newNode(project);
			left->left = getTree(project);
            left->opcode = opcode;
        }
    }   
	*lex = getLex(project);

    /* --- exponentiation */
    while ( *lex == 31 )
    {
		*lex = getLex(project);
        bracket = 0;
        if ( *lex == 1 )
        {
            bracket = 1;
			*lex = getLex(project);
        }
        if ( *lex != 7 )
        {
            project->Err = 1;
            return NULL;
        }
        right = newNode(project);
        right->opcode = *lex;
        right->fvalue = project->Fvalue;
        node = newNode(project);
        node->left = left;
        node->right = right;
        node->opcode = 31;
        left = node;
        if (bracket)
        {
			*lex = getLex(project);
            if ( *lex != 2 )
            {
                project->Err = 1;
                return NULL;
            }
        }
		*lex = getLex(project);
    }
    return left;
}

//=============================================================================

ExprTree * getOp(Project* project, int *lex)
{
    int opcode;
    ExprTree *left;
    ExprTree *right;
    ExprTree *node;
    int neg = 0;

    *lex = getLex(project);
    if (project->PrevLex == 0 || project->PrevLex == 1)
    {
        if ( *lex == 4 )
        {
            neg = 1;
			*lex = getLex(project);
        }
		else if (*lex == 3) *lex = getLex(project);
    }
	left = getSingleOp(project,lex);
    while ( *lex == 5 || *lex == 6 )
    {
        opcode = *lex;
		*lex = getLex(project);
		right = getSingleOp(project,lex);
	node = newNode(project);
	if (project->Err) return NULL;
        node->left = left;
        node->right = right;
        node->opcode = opcode;
        left = node;
    }
    if ( neg )
    {
        node = newNode(project);
        if (project->Err) return NULL;
        node->left = left;
        node->right = NULL;
        node->opcode = 9;
        left = node;
    }
    return left;
}

//=============================================================================

ExprTree * getTree(Project* project)
{
    int      lex;
    int      opcode;
    ExprTree *left;
    ExprTree *right;
    ExprTree *node;

	left = getOp(project,&lex);
    for (;;)
    {
        if ( lex == 0 || lex == 2 )
        {
            if ( lex == 2 ) project->Bc--;
            break;
        }

        if (lex != 3 && lex != 4 )
        {
            project->Err = 1;
            break;
        }

        opcode = lex;
		right = getOp(project, &lex);
	node = newNode(project);
	if (project->Err) break;
        node->left = left;
        node->right = right;
        node->opcode = opcode;
        left = node;
    } 
    return left;
}

//=============================================================================

void traverseTree(ExprTree *tree, MathExpr **expr)
// Converts binary tree to linked list (postfix format)
{
    MathExpr *node;
    if ( tree == NULL) return;
    traverseTree(tree->left,  expr);
    traverseTree(tree->right, expr);
    node = (MathExpr *) malloc(sizeof(MathExpr));
    node->fvalue = tree->fvalue;
    node->opcode = tree->opcode;
    node->ivar = tree->ivar;
    node->next = NULL;
    node->prev = (*expr);
    if (*expr) (*expr)->next = node;
    (*expr) = node;
}

//=============================================================================

void deleteTree(ExprTree *tree)
{
    if (tree)
    {
        if (tree->left)  deleteTree(tree->left);
        if (tree->right) deleteTree(tree->right);
        free(tree);
    }
}

//=============================================================================

// Turn on "precise" floating point option                                     //(5.1.008)
#pragma float_control(precise, on, push)                                       //(5.1.008)

double mathexpr_eval(Project* project, MathExpr *expr, double(*getVariableValue) (struct Project*, int))
//  Mathematica expression evaluation using a stack
{
    
// --- Note: the ExprStack array must be declared locally and not globally
//     since this function can be called recursively.

    double ExprStack[MAX_STACK_SIZE];
    MathExpr *node = expr;
    double r1, r2;
    int stackindex = 0;
    
    ExprStack[0] = 0.0;
    while(node != NULL)
    {
	switch (node->opcode)
	{
	    case 3:  
		r1 = ExprStack[stackindex];
		stackindex--;
		r2 = ExprStack[stackindex];
		ExprStack[stackindex] = r2 + r1;
		break;

        case 4:  
		r1 = ExprStack[stackindex];
		stackindex--;
		r2 = ExprStack[stackindex];
		ExprStack[stackindex] = r2 - r1;
		break;

        case 5:  
		r1 = ExprStack[stackindex];
		stackindex--;
		r2 = ExprStack[stackindex];
		ExprStack[stackindex] = r2 * r1;
		break;

        case 6:  
		r1 = ExprStack[stackindex];
		stackindex--;
		r2 = ExprStack[stackindex];
		ExprStack[stackindex] = r2 / r1;
		break;				

        case 7:  
		stackindex++;
		ExprStack[stackindex] = node->fvalue;
		break;

        case 8:
        if (getVariableValue != NULL)
        {
           r1 = getVariableValue(project, node->ivar);
        }
        else r1 = 0.0;
		stackindex++;
		ExprStack[stackindex] = r1;
		break;

        case 9: 
		ExprStack[stackindex] = -ExprStack[stackindex];
		break;

        case 10: 
		r1 = ExprStack[stackindex];
		r2 = cos(r1);
		ExprStack[stackindex] = r2;
		break;

        case 11: 
		r1 = ExprStack[stackindex];
		r2 = sin(r1);
		ExprStack[stackindex] = r2;
		break;

        case 12: 
		r1 = ExprStack[stackindex];
		r2 = tan(r1);
		ExprStack[stackindex] = r2;
		break;

        case 13: 
		r1 = ExprStack[stackindex];
		if (r1 == 0.0) r2 = 0.0;
		else r2 = 1.0/tan( r1 );    
		ExprStack[stackindex] = r2;
		break;

        case 14: 
		r1 = ExprStack[stackindex];
		r2 = fabs( r1 );       
		ExprStack[stackindex] = r2;
		break;

        case 15: 
		r1 = ExprStack[stackindex];
		if (r1 < 0.0) r2 = -1.0;
		else if (r1 > 0.0) r2 = 1.0;
		else r2 = 0.0;
		ExprStack[stackindex] = r2;
		break;

        case 16: 
		r1 = ExprStack[stackindex];
		if (r1 < 0.0) r2 = 0.0;
		else r2 = sqrt( r1 );     
		ExprStack[stackindex] = r2;
		break;

        case 17: 
		r1 = ExprStack[stackindex];
		if (r1 <= 0) r2 = 0.0;
		else r2 = log(r1);
		ExprStack[stackindex] = r2;
		break;

        case 18: 
		r1 = ExprStack[stackindex];
		r2 = exp(r1);
		ExprStack[stackindex] = r2;
		break;

        case 19: 
		r1 = ExprStack[stackindex];
		r2 = asin( r1 );
		ExprStack[stackindex] = r2;
		break;

        case 20: 
		r1 = ExprStack[stackindex];
		r2 = acos( r1 );      
		ExprStack[stackindex] = r2;
		break;

        case 21: 
		r1 = ExprStack[stackindex];
		r2 = atan( r1 );      
		ExprStack[stackindex] = r2;
		break;

        case 22: 
		r1 = ExprStack[stackindex];
		r2 = 1.57079632679489661923 - atan(r1);  
		ExprStack[stackindex] = r2;
		break;

        case 23:
		r1 = ExprStack[stackindex];
		r2 = (exp(r1)-exp(-r1))/2.0;
		ExprStack[stackindex] = r2;
		break;

        case 24: 
		r1 = ExprStack[stackindex];
		r2 = (exp(r1)+exp(-r1))/2.0;
		ExprStack[stackindex] = r2;
		break;

        case 25: 
		r1 = ExprStack[stackindex];
		r2 = (exp(r1)-exp(-r1))/(exp(r1)+exp(-r1));
		ExprStack[stackindex] = r2;
		break;

        case 26: 
		r1 = ExprStack[stackindex];
		r2 = (exp(r1)+exp(-r1))/(exp(r1)-exp(-r1));
		ExprStack[stackindex] = r2;
		break;

        case 27: 
		r1 = ExprStack[stackindex];
		if (r1 == 0.0) r2 = 0.0;
		else r2 = log10( r1 );     
		ExprStack[stackindex] = r2;
		break;

        case 28:
 		r1 = ExprStack[stackindex];
		if (r1 <= 0.0) r2 = 0.0;
		else           r2 = 1.0;
		ExprStack[stackindex] = r2;
		break;
               
        case 31: 
		r1 = ExprStack[stackindex];
		r2 = ExprStack[stackindex-1];
		if (r2 <= 0.0) r2 = 0.0;
		else r2 = exp(r1*log(r2));
		ExprStack[stackindex-1] = r2;
		stackindex--;
		break;
        }
        node = node->next;
    }
    r1 = ExprStack[stackindex];

    // Set result to 0 if it is NaN due to an illegal math op                  //(5.1.008)
    if ( r1 != r1 ) r1 = 0.0;                                                  //(5.1.008)

    return r1;
}

// Turn off "precise" floating point option                                    //(5.1.008)
#pragma float_control(pop)                                                     //(5.1.008)

//=============================================================================

void mathexpr_delete(MathExpr *expr)
{
    if (expr) mathexpr_delete(expr->next);
    free(expr);
}

//=============================================================================

MathExpr * mathexpr_create(Project*  project, char *formula, int(*getVar) (struct Project*, char *))
{
    ExprTree *tree;
    MathExpr *expr = NULL;
    MathExpr *result = NULL;
    getVariableIndex = getVar;
    project->Err = 0;
    project->PrevLex = 0;
    project->CurLex = 0;
    project->S = formula;
    project->Len = strlen(project->S);
     project->Pos = 0;
    project->Bc = 0;
	tree = getTree(project);
    if (project->Bc == 0 && project->Err == 0)
    {
	    traverseTree(tree, &expr);
	    while (expr)
	    {
            result = expr;
            expr = expr->prev;
        }
    }
    deleteTree(tree);
    return result;
}
