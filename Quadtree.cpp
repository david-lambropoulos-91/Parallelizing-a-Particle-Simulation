#include<stdio.h>
#include<stdlib.h>

//typedef struct Point Point;
struct Point
{
int x;
int y;
};

//typedef struct Node Node;
struct Node
{
int posX;
int posY;
int width;
int height;
Node *child[4];         //Changed to Node *child[4] rather than Node ** child[4]
Point pointArray[5000];
};

void BuildQuadTree(Node *n);
void PrintQuadTree(Node *n, int depth);
void DeleteQuadTree(Node *n);
Node *BuildNode(Node *n, Node  *nParent, int index);

void setnode(Node *xy, int x, int y, int w, int h)
{
int i;
xy->posX = x;
xy->posY = y;
xy->width = w;
xy->height = h;
for (i = 0; i < 5000; i++)
{
xy->pointArray[i].x = 560;
xy->pointArray[i].y = 560;
}
//Initialises child-nodes to NULL - better safe than sorry
for (int i = 0; i < 4; i++)
xy->child[i] = NULL;
}

int randn()
{
int a;
a = rand() % 501;
return a;
}

int pointArray_size(Node *n)
{
int m = 0, i;
for (i = 0; i <= 5000; i++)
if (n->pointArray[i].x <= 500 && n->pointArray[i].y <= 500)
m++;
return (m + 1);
}

//int main()
//{
// Initialize the root node
//Node * rootNode = new Node;    //Initialised node
//int i, x[5000], y[5000];
//setnode(rootNode, 0, 0, 500, 500);
//for (i = 0; i < 5000; i++)
//{
//x[i] = randn();
//y[i] = randn();
//}
//for (i = 0; i < 5000; i++)
//{
//rootNode->pointArray[i].x = x[i];
//rootNode->pointArray[i].y = y[i];
//}
// Create the quadTree
//BuildQuadTree(rootNode);
//PrintQuadTree(rootNode, 0);    //Added function to print for easier debugging
//DeleteQuadTree(rootNode);
//return 0;
//}

void BuildQuadTree(Node *n)
{
Node * nodeIn = new Node;   //Initialised node
int points = pointArray_size(n);
if (points > 100)
{
for (int k = 0; k < 4; k++)
{
n->child[k] = new Node;      //Initialised node
nodeIn = BuildNode(n->child[k], n, k);
BuildQuadTree(nodeIn);
}
}
}

void PrintQuadTree(Node *n, int depth)
{
for (int i = 0; i < depth; i++)
printf("\t");
if (n->child[0] == NULL)
{
int points = pointArray_size(n);
printf("Points: %d\n", points);
return;
}
else if (n->child[0] != NULL)
{
printf("Children:\n");
for (int i = 0; i < 4; i++)
PrintQuadTree(n->child[i], depth + 1);
return;
}
}

void DeleteQuadTree(Node *n)
{
if (n->child[0] == NULL)
{
delete n;
return;
}
else if (n->child[0] != NULL)
{
for (int i = 0; i < 4; i++)
DeleteQuadTree(n->child[i]);
return;
}
}

Node *BuildNode(Node *n, Node *nParent, int index)
{
int numParentPoints, i, j = 0;
// 1) Creates the bounding box for the node
// 2) Determines which points lie within the box
/*
Position of the child node, based on index (0-3), is determined in this order:
| 1 | 0 |
| 2 | 3 |*/

setnode(n, 0, 0, 0, 0);
switch (index)
{
case 0: // NE
n->posX = nParent->posX + nParent->width / 2;
n->posY = nParent->posY + nParent->height / 2;
break;

case 1: // NW
n->posX = nParent->posX;
n->posY = nParent->posY + nParent->height / 2;
break;

case 2: // SW
n->posX = nParent->posX;
n->posY = nParent->posY;
break;

case 3: // SE

n->posX = nParent->posX + nParent->width / 2;
n->posY = nParent->posY;
break;
}
// Width and height of the child node is simply 1/2 of the parent node's width and height
n->width = nParent->width / 2;
n->height = nParent->height / 2;
// The points within the child node are also based on the index, similiarily to the position
numParentPoints = pointArray_size(nParent);
switch (index)
{
case 0: // NE
for (i = 0; i < numParentPoints - 1; i++)
{
// Check all parent points and determine if it is in the top right quadrant
if (nParent->pointArray[i].x <= 500 && nParent->pointArray[i].x > nParent->posX + nParent->width / 2 && nParent->pointArray[i].y > nParent->posY + nParent->height / 2 && nParent->pointArray[i].x <= nParent->posX + nParent->width && nParent->pointArray[i].y <= nParent->posY + nParent->height)
{
// Add the point to the child node's point array
n->pointArray[j].x = nParent->pointArray[i].x;
n->pointArray[j].y = nParent->pointArray[i].y;
j++;
}
}
break;
case 1: // NW
for (i = 0; i < numParentPoints - 1; i++)
{
// Check all parent points and determine if it is in the top left quadrant
if (nParent->pointArray[i].x <= 500 && nParent->pointArray[i].x > nParent->posX && nParent->pointArray[i].y > nParent->posY + nParent->height / 2 && nParent->pointArray[i].x <= nParent->posX + nParent->width / 2 && nParent->pointArray[i].y <= nParent->posY + nParent->height)
{
// Add the point to the child node's point array
n->pointArray[j].x = nParent->pointArray[i].x;
n->pointArray[j].y = nParent->pointArray[i].y;
j++;
}
}
break;
case 2: // SW
for (i = 0; i < numParentPoints - 1; i++)
{
// Check all parent points and determine if it is in the bottom left quadrant
if (nParent->pointArray[i].x <= 500 && nParent->pointArray[i].x > nParent->posX && nParent->pointArray[i].y > nParent->posY && nParent->pointArray[i].x <= nParent->posX + nParent->width / 2 && nParent->pointArray[i].y <= nParent->posY + nParent->height / 2)
{
// Add the point to the child node's point array
n->pointArray[j].x = nParent->pointArray[i].x;
n->pointArray[j].y = nParent->pointArray[i].y;
j++;
}
}
break;
case 3: // SE
for (i = 0; i < numParentPoints - 1; i++)
{
// Check all parent points and determine if it is in the bottom right quadrant
if (nParent->pointArray[i].x <= 500 && nParent->pointArray[i].x > nParent->posX + nParent->width / 2 && nParent->pointArray[i].y > nParent->posY && nParent->pointArray[i].x <= nParent->posX + nParent->width && nParent->pointArray[i].y <= nParent->posY + nParent->height / 2)
{
// Add the point to the child node's point array
n->pointArray[j].x = nParent->pointArray[i].x;
n->pointArray[j].y = nParent->pointArray[i].y;
j++;
}
}
break;
}
return n;
}
