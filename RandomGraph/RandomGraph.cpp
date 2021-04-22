// RandomGraph.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include<iostream>
#include <unordered_set>
#include <list>
#include <queue>
#include <time.h>
#include <iostream>       
#include <thread> 
#include <iostream>
#include <fstream>

//Netanel Mayo 
using namespace std;

struct Vertice
{
	int data;
	int deg = 0;
	struct Vertice* next;
};

Vertice** Graph;
Vertice* new_Vertice;

int E = 0;
const int V = 1000;
double P = 0.5;


Vertice** build_random_graph(int v, double p);
void display(Vertice** Graph, int V);
void insert(Vertice** Graph, int v_data);
int Connectivity(Vertice** graph);
int Is_Isolated(Vertice** graph);
int BFS(int init, Vertice** graph, int n, int diameter[]);
int FindDiameter(Vertice** graph);
void P_arrayCalculate(double median, double(&returnArray)[10]);
void GenerateGraph(double P, int(&count_Isolated), int(&count_Connectivity), int(&count_Diameter), string p_test);
void FreePointers(Vertice** Graph);
void SaveResultToFile(double P_threshold[], double P_threshold_result[], string p_test);
void GenerateLoopsGraph(double P_threshold[], double P_threshold_result[], int amountOfGraphGenerate, string p_test);


int main()
{
	srand((unsigned)time(0));
	
	int count_Isolated = 0;
	int count_Connectivity = 0;
	int count_Diameter = 0;

	double P_threshold1_connectivity[10] = { 0 };
	double P_threshold1_connectivity_result[10] = { 0 };

	double P_threshold2_diameter[10] = { 0 };
	double P_threshold2_diameter_result[10] = { 0 };

	double P_threshold3_isolated[10] = { 0 };
	double P_threshold3_isolated_result[10] = { 0 };

	double threshold1_connectivity = (log(V) / V);
	P_arrayCalculate(threshold1_connectivity, P_threshold1_connectivity);

	double threshold2_diameter = sqrt((2 * log(V)) / (V));
	P_arrayCalculate(threshold2_diameter, P_threshold2_diameter);

	double threshold3_isolated = threshold1_connectivity;
	P_arrayCalculate(threshold3_isolated, P_threshold3_isolated);

	int amountOfGraphGenerate =500;
	GenerateLoopsGraph(P_threshold1_connectivity, P_threshold1_connectivity_result, amountOfGraphGenerate, "connectivity");
	GenerateLoopsGraph(P_threshold2_diameter, P_threshold2_diameter_result, amountOfGraphGenerate, "diameter");
	GenerateLoopsGraph(P_threshold3_isolated, P_threshold3_isolated_result, amountOfGraphGenerate, "isolated");

	return 0;
}

//This Func if for call each test that build 500 times random graph for each P  
void GenerateLoopsGraph(double P_threshold[], double P_threshold_result[], int amountOfGraphGenerate, string p_test) {

	int count_Isolated = 0;
	int count_Connectivity = 0;
	int count_Diameter = 0;

	// current date/time based on current system
	time_t now = time(0);
	// convert now to string form
	char* dt = ctime(&now);

	std::cout << "\n"<< dt <<"Start Process P test: " << p_test<<endl;
	for (int p = 0; p < 10; p++)
	{
		for (int i = 0; i < amountOfGraphGenerate; i++)
		{
			GenerateGraph(P_threshold[p], count_Isolated, count_Connectivity, count_Diameter, p_test);	
		}
		if (p_test == "isolated") {
			double sum = ((double)count_Isolated) / (double)amountOfGraphGenerate;
			P_threshold_result[p] = sum;
		}
		else if (p_test == "connectivity") {
			double sum = ((double)count_Connectivity) / (double)amountOfGraphGenerate;
			P_threshold_result[p] = (double)sum;
		}
		else if (p_test == "diameter") {
			double sum = ((double)count_Diameter) / (double)amountOfGraphGenerate;;
			P_threshold_result[p] = sum;}
		
		count_Isolated = count_Connectivity = count_Diameter = 0;
	}
	for (int i = 0; i < 10; i++)
	{
		printf(" , %f", P_threshold[i]);
	}
	std::cout << endl;

	for (int i = 0; i < 10; i++)
	{
		printf(" , %f", P_threshold_result[i]);
	}

	now = time(0);
	dt = ctime(&now);

	std::cout << "\n"<< dt <<"End Process P test: " << p_test <<endl;
	try {
		
		SaveResultToFile(P_threshold, P_threshold_result, p_test);
	}
	catch (exception mess) {
		std::cout << "\nError To Save file results P test: " << p_test << endl;
	}	
}

//this func is for saveing the results from each test
void SaveResultToFile(double P_threshold[], double P_threshold_result[], string p_test)
{

	// file pointer
	fstream myfile;
	// opens an existing csv file or creates a new file.
	myfile.open("results.csv", ios::out | ios::app);
	myfile << "P " << p_test << ", ";
	for (int i = 0; i < 10; i++)
	{
		myfile << P_threshold[i] << ", ";
	}
	myfile << "\n";
	myfile << "Probability that the feature exists" << ", ";
	for (int i = 0; i < 10; i++)
	{
		myfile << P_threshold_result[i] << ", ";
	}
	myfile <<"\n";
	myfile.close();
}

//this func is for calling to func that generate random graph and for calling to the func that  calculate isolated, connectivity and daimeter 
void GenerateGraph(double P, int(&count_Isolated), int(&count_Connectivity), int(&count_Diameter), string p_test)
{
	Vertice** Graph;
	Graph = build_random_graph(V, P);
	//display(Graph, V);
	if (Is_Isolated(Graph) == 0)
	{
		//printf("Graph Isolated \n");
		count_Isolated++;
		if (p_test=="connectivity" && Connectivity(Graph) == 1)
		{			
			//printf("Graph Connectivity \n");
			count_Connectivity++;
		}
		else {
			//printf("Graph Not Connectivity \n");
		}
		if (p_test == "diameter" && FindDiameter(Graph) <= 2)
		{		
			count_Diameter++;
		}
		//printf("Diameter of Graph is %d\n", FindDiameter(Graph));
	}
	else {
		/*printf("Graph Not Isolated \n");
		printf("Graph Not Connectivity \n");
		printf("Diameter of Graph is %d\n", INT_MAX);*/
	}
	
	//printf("Is_Isolated Time stamp %d,\n", after - before);
	FreePointers(Graph);
}

//This func is for calculate dynamic P array for referance
void P_arrayCalculate(double median, double(&returnArray)[10])
{
	//double num = (median / 6);

	double num1 = (median / 6);
	//double num2 = ((1 - median) / 6);
	double num2 = ((1 - median) / 60);
	//printf("median %f", median);
	returnArray[0] = num1;
	//printf(", %f", returnArray[0]);
	for (int i = 1; i < 10; i++)
	{
		if (i >= 5)
		{
			//num2 = num2 + returnArray[i - 1];
			returnArray[i] = (num2 + returnArray[i - 1]);
		}
		else {
			returnArray[i] = (num1 + returnArray[i - 1]);
		}
		//printf(", %f", returnArray[i]);
	}
}

//This func is for build random graph
Vertice** build_random_graph(int v, double p)
{
	E = 0;
	Vertice** Graph = new Vertice * [v];

	for (int i = 0; i < v; i++)
	{
		new_Vertice = new Vertice();
		new_Vertice->data = i;
		new_Vertice->deg = 0;
		new_Vertice->next = NULL;
		Graph[i] = new_Vertice;
	}

	double probability;
	//double probability2;
	int count = 0;
	for (int i = 0; i < v; i++)
	{
		for (int y = i + 1; y < v; y++)
		{
			probability = (rand() / (double)(RAND_MAX + 1));
			if (probability <= p)
			{
				insert(&Graph[i], y);
				insert(&Graph[y], i);
				E++;
			}
		}
	}
	return Graph;
}

//this func is for insert new connection between 2 Vertices 
void insert(Vertice** Graph, int v_data)
{
	struct Vertice* ptr;
	struct Vertice* temp;
	new_Vertice = new Vertice();
	new_Vertice->data = v_data;
	new_Vertice->next = NULL;
	new_Vertice->deg = 0;


	//int size = (*Graph)->deg;

	(*Graph)->deg++;
	ptr = *Graph;

	if ((ptr)->next == NULL) {
		ptr->next = new_Vertice;
	}
	else {
		temp=ptr->next;
		new_Vertice->next = temp;
		ptr->next = new_Vertice;
		/*while (ptr->next != NULL)
			ptr = ptr->next;
		ptr->next = new_Vertice;*/
	}
}

//for debug only - this is show to the graph 
void display(Vertice** Graph, int V)
{
	struct Vertice* ptr = new Vertice;
	for (int i = 0; i < V; i++)
	{
		ptr = Graph[i];
		while (ptr != NULL)
		{
			cout << ptr->data << ", ";
			ptr = ptr->next;
		}
		printf("\n");
	}
}

//Frees memory for pointers
void FreePointers(Vertice** Graph)
{
	Vertice* temp;
	Vertice** temp_Graph;
	for (int i = 0; i < V; i++)
	{
		while (Graph[i] != NULL)
		{
			temp = Graph[i];
			Graph[i] = Graph[i]->next;
			delete temp;
		}
	}
	if (Graph != NULL) {
		delete[] Graph;
	}
}

//Calculate Connectivity if it found return 1 else return 0, use BFS algorithm 
int Connectivity(Vertice** graph)
{
	// Initializing queue
	struct Vertice* ptr;
	queue<int> q;
	int visited[V] = { 0 };

	int init = 0;
	// Pushing each node in queue
	q.push(init);

	// Mark the traversed node visited
	visited[init] = 1;
	while (!q.empty()) {
		int u = q.front();
		ptr = graph[u];
		q.pop();
		while (ptr != NULL)
		{
			int data = ptr->data;
			if (visited[data] == 0) {
				visited[data] = 1;
				q.push(data);
			}
			ptr = ptr->next;
		}
	}

	for (int i = 0; i < V; i++)
	{
		if (visited[i] == 0) {
			return 0;
		}
	}
	return 1;
}

//Checks if there is at least one single node if yes->return 1 else return 0
int Is_Isolated(Vertice** graph)
{
	for (int i = 0; i < V; i++) {
		if (graph[i]->deg == 0)
		{
			return 1;
		}
	}
	return 0;
}

//Calculate and returns the diameter of the graph use BFS algorithm 
int FindDiameter(Vertice** graph)
{
	int* diameter = new int[V];
	int init = BFS(0, graph, V, diameter);
	int val = BFS(init, graph, V, diameter);	
	int diam = diameter[val];
	delete[] diameter;
	return diam;
}

//BFS algorithm 
int BFS(int init, Vertice** graph, int n, int diameter[])
{
	// Initializing queue
	struct Vertice* ptr;
	//ptr = graph[u];

	queue<int> q;
	//q.push(init);

	int visited[V] = { 0 };

	for (int i = 0; i < n; i++) {
		//visited[i] = 0;
		diameter[i] = 0;
	}

	// Pushing each node in queue
	q.push(init);

	// Mark the traversed node visited
	visited[init] = 1;
	while (!q.empty()) {
		int u = q.front();
		ptr = graph[u];
		q.pop();

		while (ptr != NULL)
		{
			int data = ptr->data;

			if (visited[data] == 0) {
				visited[data] = 1;

				// Considering weight of edges equal to 1
				diameter[data] += diameter[u] + 1;
				q.push(data);
			}
			ptr = ptr->next;
		}
	}
	// return index of max value in diameter
	return int(max_element(diameter + 1,
		diameter + n + 1)
		- diameter);
}

