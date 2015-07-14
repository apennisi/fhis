/*
 *  FHIS Image Segmentation Library
 *
 *
 *  Written by Andrea Pennisi
 *
 *  Please, report suggestions/comments/bugs to
 *  andrea.pennisi@gmail.com
 *
 */

#ifndef MERGETRIANGLES_H
#define MERGETRIANGLES_H

#include <math.h>
#include <opencv2/opencv.hpp>
#include <utils.h>

class MergeTriangles
{
    typedef struct edge_tag {

         int idx1, idx2;

         double merge_cost;
         double total_length;

         struct edge_tag *next1, *prev1;
         struct edge_tag *next2, *prev2;

         int heap_index;
     } edge;

     typedef struct {
         struct edge_tag *edge_list;
     } node;

     struct {
         int count;
         node *buf;
     } node_buffer;

     struct {
         int count;
         edge *buf;
     } edge_buffer;

     struct {
         int count;
         edge **elts;
     } edge_heap;

public:
    MergeTriangles(const std::vector<double> &neighbors,
                   const std::vector<double> &color_weigths,
                   const std::vector<double> &edge_lengths,
                   const int &ntri);
    void getNewLabels(std::vector<double> &merge_cost_, std::vector<int> &old_labels_,
                      std::vector<int> &new_labels_, std::vector<double> &merge_lengths_);
private:
    edge **adjacency_table;

    int other_idx(edge *theEdge, const int &idx);
    void update_index (edge *theEdge, const int &old_idx, const int &new_idx);
    struct edge_tag *get_next_edge (edge *theEdge, const int &idx);
    struct edge_tag *get_prev_edge (edge *theEdge, const int &idx);
    void set_next_edge (edge *theEdge, const int &idx, struct edge_tag *ptr);
    void set_prev_edge (edge *theEdge, const int &idx, struct edge_tag *ptr);
    void add_to_adjacency_list(const int &node_index, edge *theEdge);
    void remove_from_adjacency_list(const int &node_index, edge *theEdge);
    void heap_swap(const int &idx1, const int &idx2);
    void heapify_up(int k);
    void heapify_down(int k);
    void add_to_edge_heap(edge *theEdge);
    void remove_from_edge_heap(edge *theEdge);
    void InitEdges(const std::vector<double> &SN, const std::vector<double> &edge_weights,
                   const std::vector<double> &edge_lengths);
    void merge_edges (edge *edge1, edge *edge2);
    void merge_nodes(const int &idx1, const int &idx2);
    void InitNodes();
    int sizeL;
    int ntri;
    int nedges;
    std::vector<double> merge_cost, merge_length;
    std::vector<int> old_labels, new_labels;
};

#endif // MERGETRIANGLES_H
