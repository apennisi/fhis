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

#include "mergetriangles.h"

MergeTriangles::MergeTriangles(const std::vector<double> &neighbors,
                               const std::vector<double> &color_weigths,
                               const std::vector<double> &edge_lengths,
                               const int &ntri) {
    int i;
    edge *best_edge;

    this->ntri = ntri;

    nedges = ((ntri * 3) / 2) + 1;

    node_buffer.count = 0;
    node_buffer.buf = new node[ntri];

    edge_buffer.count = 0;
    edge_buffer.buf = new edge[nedges];

    edge_heap.count = 0;
    edge_heap.elts = new edge*[nedges];

    adjacency_table = new edge*[ntri];

    for(int i = 0; i < ntri; i++) {
        *(adjacency_table + i) = 0;
    }

    // Initialize the nodes
    InitNodes();

    // Initialize the edges
    InitEdges(neighbors, color_weigths, edge_lengths);

    sizeL = ntri-1;

    new_labels.resize(sizeL);
    old_labels.resize(sizeL);
    merge_cost.resize(sizeL);
    merge_length.resize(sizeL);

    for (i=0; i < sizeL; ++i) {

        best_edge = edge_heap.elts[0];

        new_labels[i]   = best_edge->idx1 + 1.;
        old_labels[i]   = best_edge->idx2 + 1.;
        merge_cost[i]   = best_edge->merge_cost;
        merge_length[i] = best_edge->total_length;

        remove_from_adjacency_list (best_edge->idx1, best_edge);
        remove_from_adjacency_list (best_edge->idx2, best_edge);

        remove_from_edge_heap(best_edge);

        merge_nodes(best_edge->idx1, best_edge->idx2);
    }

    delete[] node_buffer.buf;
    delete[] edge_buffer.buf;
    delete[] edge_heap.elts;
    delete[] adjacency_table;
}

void MergeTriangles::InitNodes() {
    node *node_ptr;
    int i;

    for (i=0, node_ptr = node_buffer.buf; i < ntri; ++i, ++node_ptr) {
        node_ptr->edge_list = NULL;
    }
}

int MergeTriangles::other_idx(edge *theEdge, const int &idx) {
    return ( (theEdge->idx1 == idx) ? theEdge->idx2 : theEdge->idx1 );
}

void MergeTriangles::update_index (edge *theEdge, const int &old_idx, const int &new_idx) {
    if (theEdge->idx1 == old_idx)
        theEdge->idx1 = new_idx;
    else
        theEdge->idx2 = new_idx;
}

struct MergeTriangles::edge_tag *MergeTriangles::get_next_edge (edge *theEdge, const int &idx) {
    return (theEdge->idx1 == idx) ? theEdge->next1 : theEdge->next2;
}

struct MergeTriangles::edge_tag *MergeTriangles::get_prev_edge (edge *theEdge, const int &idx) {
    return (theEdge->idx1 == idx) ? theEdge->prev1 : theEdge->prev2;
}


void MergeTriangles::set_next_edge (edge *theEdge, const int &idx, struct edge_tag *ptr) {
    if (theEdge->idx1 == idx) {
        theEdge->next1 = ptr;
    } else {
        theEdge->next2 = ptr;
    }
}

void MergeTriangles::set_prev_edge (edge *theEdge, const int &idx, struct edge_tag *ptr) {
    if (theEdge->idx1 == idx) {
        theEdge->prev1 = ptr;
    } else {
        theEdge->prev2 = ptr;
    }
}

void MergeTriangles::add_to_adjacency_list(const int &node_index, edge *theEdge) {

    node *theNode = node_buffer.buf + node_index;
    edge *head_elt = theNode->edge_list;

    set_next_edge (theEdge, node_index, head_elt);
    set_prev_edge (theEdge, node_index, NULL);

    if (head_elt) {
        set_prev_edge (head_elt, node_index, theEdge);
    }

    theNode->edge_list = theEdge;
}

void MergeTriangles::remove_from_adjacency_list(const int &node_index, edge *theEdge) {
    node *theNode = node_buffer.buf + node_index;
    edge *prev = get_prev_edge (theEdge, node_index);
    edge *next = get_next_edge (theEdge, node_index);

    if (prev) {
        set_next_edge(prev, node_index, next);
    } else {
        theNode->edge_list = next;
    }

    if (next) {
        set_prev_edge (next, node_index, prev);
    }
}


void MergeTriangles::heap_swap(const int &idx1, const int &idx2) {

    edge* temp = edge_heap.elts[idx1];
    edge_heap.elts[idx1] = edge_heap.elts[idx2];
    edge_heap.elts[idx2] = temp;

    edge_heap.elts[idx1]->heap_index = idx1;
    edge_heap.elts[idx2]->heap_index = idx2;
}


void MergeTriangles::heapify_up(int k) {

    long double merge_cost = edge_heap.elts[k]->merge_cost;
    int parent;

    while (k > 0) {
        parent = (k-1) / 2;
        if (edge_heap.elts[parent]->merge_cost > merge_cost) {
            heap_swap(parent, k);
            k = parent;
        } else
            break;
    }
}

void MergeTriangles::heapify_down(int k) {
    long double merge_cost = edge_heap.elts[k]->merge_cost;
    int child;

    while (1) {

        child = 2*k + 1;

        if (child >= edge_heap.count)
                    break;

        if ( ((child+1) < edge_heap.count) && (edge_heap.elts[child+1]->merge_cost < edge_heap.elts[child]->merge_cost) )
            child = child + 1;

        if (edge_heap.elts[child]->merge_cost < merge_cost) {
            heap_swap(child, k);
            k = child;
        } else
            break;
    }
}

void MergeTriangles::add_to_edge_heap(edge *theEdge) {
    edge_heap.elts[edge_heap.count] = theEdge;
    theEdge->heap_index = edge_heap.count;
    heapify_up(edge_heap.count);
    ++edge_heap.count;
}

void MergeTriangles::remove_from_edge_heap(edge *theEdge) {

    int idx = theEdge->heap_index;
    long double new_merge_cost, old_merge_cost;

    --(edge_heap.count);

    old_merge_cost = edge_heap.elts[idx]->merge_cost;
    new_merge_cost = edge_heap.elts[edge_heap.count]->merge_cost;

    heap_swap(idx, edge_heap.count);

    if (new_merge_cost < old_merge_cost)
        heapify_up(idx);
    else
        heapify_down(idx);
}

void MergeTriangles::InitEdges(const std::vector<double> &SN, const std::vector<double> &edge_weights,
                               const std::vector<double> &edge_lengths) {
    int i, j, idx;
    edge *theEdge;

    int offset = 0;

    for (i=0; i < ntri; ++i) {
        offset = i*3;
        for (j=0; j < 3; ++j) {
            idx = SN[offset + j];

            if (idx > 0) {

                if (idx > ntri) {
                    std::cerr << "index entry in SN too large" << std::endl;
                    exit(-1);
                }

                idx = idx - 1;

                if (i < idx) {

                    theEdge = edge_buffer.buf + edge_buffer.count;
                    ++edge_buffer.count;

                    theEdge->idx1 = i;
                    theEdge->idx2 = idx;

                    theEdge->merge_cost   = edge_weights[offset + j];
                    theEdge->total_length = edge_weights[offset + j];


                    if (theEdge->merge_cost < 0) {
                        std::cerr << "edge_weight entry negative" << std::endl;
                        exit(-1);
                    }

                    if (theEdge->total_length < 0) {
                        std::cerr << "edge_length entry is < 0" << std::endl;
                        exit(-1);
                    }

                    theEdge->prev1 = NULL;
                    theEdge->next1 = NULL;

                    theEdge->prev2 = NULL;
                    theEdge->next2 = NULL;


                    add_to_adjacency_list(i, theEdge);
                    add_to_adjacency_list(idx, theEdge);
                    add_to_edge_heap(theEdge);
                }

            }
        }
    }
}


void MergeTriangles::merge_edges (edge *edge1, edge *edge2) {
    long double merge_cost1, merge_cost2;
    long double length1, length2;
    long double new_merge_cost, new_length;

    merge_cost1 = edge1->merge_cost;
    merge_cost2 = edge2->merge_cost;

    length1 = edge1->total_length;
    length2 = edge2->total_length;

    new_length = length1 + length2;

    new_merge_cost = (long double)(length1*merge_cost1 + length2*merge_cost2) / (long double)new_length;

    edge1->merge_cost = new_merge_cost;
    edge1->total_length = new_length;

    if (new_merge_cost < merge_cost1)
        heapify_up (edge1->heap_index);
    else
        heapify_down (edge1->heap_index);
}


void MergeTriangles::merge_nodes(const int &idx1, const int &idx2) {
    int head_idx;
    node *node1, *node2;
    edge *list1, *list2, *head_elt, *clear_list;

    if (idx1 == idx2) return;

    node1 = node_buffer.buf + idx1;
    node2 = node_buffer.buf + idx2;

    list1 = node1->edge_list;
    list2 = node2->edge_list;

    node1->edge_list = NULL;

    while (list2) {
        head_elt = list2;
        list2 = get_next_edge(list2, idx2);

        head_idx = other_idx(head_elt, idx2);

        update_index (head_elt, idx2, idx1);

        adjacency_table[head_idx] = head_elt;

        add_to_adjacency_list (idx1, head_elt);

    }

    clear_list = node1->edge_list;

    while (list1) {
        head_elt = list1;
        list1 = get_next_edge(list1, idx1); // Pop list1

        head_idx = other_idx(head_elt, idx1);

        if (adjacency_table[head_idx]) {

            merge_edges (adjacency_table[head_idx], head_elt);

            remove_from_adjacency_list(head_idx, head_elt);

            remove_from_edge_heap(head_elt);

        } else {
            add_to_adjacency_list (idx1, head_elt);
        }

    }

    while (clear_list) {
        adjacency_table[other_idx(clear_list, idx1)] = NULL;
        clear_list = get_next_edge(clear_list, idx1);
    }
}

void MergeTriangles::getNewLabels(std::vector<double> &merge_cost_, std::vector<int> &old_labels_,
                                  std::vector<int> &new_labels_, std::vector<double> &merge_lengths_) {

    merge_cost_ = merge_cost;
    merge_lengths_ = merge_length;
    old_labels_ = old_labels;
    new_labels_ = new_labels;
}
