//
//  main.c
//  DNAGraphics
//
//  Created by Andrew Prosikhin on 10/12/15.
//  Copyright (c) 2015 AVP. All rights reserved.
//

#include <stdio.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lodepng.h"

#define length(x,y) (sqrt((x*x)+(y*y)))

#define IMAGE_WIDTH 200
#define IMAGE_HEIGHT 200

const char *filename = "TestTriangle.png";

typedef struct matrix {
    float a;
    float b;
    float c;
    float d;
    float e;
    float f;
    float g;
    float h;
    float i;
} Matrix;

typedef struct vector {
    float x;
    float y;
    float z;
} Vector;

typedef struct color {
    float red;
    float green;
    float blue;
} color;

typedef struct vertex {
    Vector position;
    color clr;
} vertex;

typedef struct triangle {
    vertex vertex_a;
    vertex vertex_b;
    vertex vertex_c;
} triangle;

typedef struct fragment {
    float z;
    int power;
    float color;
} Fragment;

typedef struct fragment_container {
    Fragment top_fragment;
    Fragment *second_frag;
} Fragment_Container;

typedef struct fragment_buffer {
    Fragment_Container fragment_array[IMAGE_WIDTH][IMAGE_HEIGHT];
} fragment_buffer;



Matrix * new_unit() {
    Matrix *m = (Matrix *)malloc(sizeof(Matrix));
    m->a = 1.0;
    m->e = 1.0;
    m->i = 1.0;
    return m;
}

Vector * affine_transformation(Matrix *m, Vector *v) {
    Vector *r = (Vector *)malloc(sizeof(Vector));
    r->x = m->a*v->x + m->d*v->y + m->g*v->z;
    r->y = m->b*v->x + m->e*v->y + m->h*v->z;
    r->z = v->z;
    return r;
}

Vector * coord_to_bar(Matrix *m, Vector *v) {
    Vector *r = (Vector *)malloc(sizeof(Vector));
    r->x = m->a*v->x + m->d*v->y + m->g*v->z;
    r->y = m->b*v->x + m->e*v->y + m->h*v->z;
    r->z = 1 - r->x - r->y;
    return r;
}

Matrix * matrix_mult_matrix(Matrix *l, Matrix *r) {
    Matrix *m = (Matrix *)malloc(sizeof(Matrix));
    m->a = l->a*r->a + l->d*r->b + l->g*r->c;
    m->b = l->b*r->a + l->e*r->b + l->h*r->c;
    m->c = l->c*r->a + l->f*r->b + l->i*r->c;
    m->d = l->a*r->d + l->d*r->e + l->g*r->f;
    m->e = l->b*r->d + l->e*r->e + l->h*r->f;
    m->f = l->c*r->d + l->f*r->e + l->i*r->f;
    m->g = l->a*r->g + l->d*r->h + l->g*r->i;
    m->h = l->b*r->g + l->e*r->h + l->h*r->i;
    m->i = l->c*r->g + l->f*r->h + l->i*r->i;
    return m;
}

void print_matrix(Matrix *m) {
    printf("%f %f %f\n%f %f %f\n%f %f %f\n", m->a, m->d, m->g, m->b, m->e, m->h, m->c, m->f, m->i);
}

Matrix * compute_converter_barycentric(Vector *a, Vector *b, Vector *c) {
    Matrix *m = (Matrix *)malloc(sizeof(Matrix));
    Matrix *t = new_unit();
    t->c = 1.0;
    t->f = 1.0;
    
    
    float d = 1/(a->x*(b->y-c->y) - b->x*(a->y-c->y) + c->x*(a->y-b->y));
    
    m->a = d*(b->y - c->y);
    m->d = d*(c->x - b->x);
    m->g = d*(b->x*c->y - c->x*b->y);
    m->b = d*(c->y - a->y);
    m->e = d*(a->x - c->x);
    m->f = d*(b->x - a->x);
    m->c = d*(a->y-b->y);
    m->h = d*(c->x*a->y - a->x*c->y);
    m->i = d*(a->x*b->y - b->x*a->y);
    return matrix_mult_matrix(t, m);
}

void printf_vector(Vector *v) {
    printf("%f\n%f\n%f\n", v->x, v->y, v->z);
}

int max(int a, int b) {
    if (a > b) {
        return a;
    }
    return b;
}

int min(int a, int b) {
    if (a > b) {
        return b;
    }
    return a;
}

float float_min(float a, float b, float c) {
    if (a > b) {
        if (b > c) {
            return c;
        } else {
            return b;
        }
    } else {
        if (a > c) {
            return c;
        } else {
            return b;
        }
    }
}

float float_max(float a, float b, float c) {
    if (a > b) {
        if (a > c) {
            return a;
        } else {
            return c;
        }
    } else {
        if (b > c) {
            return b;
        } else {
            return c;
        }
    }
}

int process_triangle(fragment_buffer *buffer, triangle *triangle) {
    int min_x = (int)float_min(triangle->vertex_a.position.x, triangle->vertex_b.position.x, triangle->vertex_c.position.x);
    int min_y = (int)float_min(triangle->vertex_a.position.y, triangle->vertex_b.position.y, triangle->vertex_c.position.y);
    int max_x = (int)float_max(triangle->vertex_a.position.x, triangle->vertex_b.position.x, triangle->vertex_c.position.x) + 1;
    int max_y = (int)float_max(triangle->vertex_a.position.y, triangle->vertex_b.position.y, triangle->vertex_c.position.y) + 1;

    if (min(min_x, min_y) < 0 || max_x > IMAGE_HEIGHT || max_y > IMAGE_HEIGHT) {
        printf("Triangle out of range!");
        return 0;
    }
    
    Matrix *m = compute_converter_barycentric(&triangle->vertex_a.position, &triangle->vertex_b.position, &triangle->vertex_c.position);
    
    Vector point, *point_baryocentric;
    int pixel_power;
    float pixel_z = 0;
    point.z = 1.0;
    for(int x = 0; x < IMAGE_WIDTH; x++)
    {
        for(int y = 0; y < IMAGE_HEIGHT; y++)
        {
            
        }
    }
    
    return 1;
}

fragment_buffer * new_fragment_buffer() {
    fragment_buffer *new_buffer = (fragment_buffer *)malloc(sizeof(fragment_buffer) + sizeof(Fragment_Container)*(IMAGE_WIDTH*IMAGE_HEIGHT - 1));
    return new_buffer;
}

void add_fragment_to_buffer(fragment_buffer *buffer, int x, int y, float z, int power) {
    Fragment_Container frag_cont;
    frag_cont.top_fragment.z = z;
    frag_cont.top_fragment.power = power;
    //frag.color; //Add Color!
    buffer->fragment_array[x][y] = frag_cont; //Primitive temporary solution!
}

void encodeOneStep(const char* filename, const unsigned char* image, unsigned width, unsigned height)
{
    /*Encode the image*/
    unsigned int error = lodepng_encode32_file(filename, image, width, height);
    
    /*if there's an error, display it*/
    if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
}

int main(int argc, const char * argv[]) {
    fragment_buffer *f_buffer = new_fragment_buffer();
    
    
    Vector *a = (Vector *)malloc(sizeof(Vector));
    Vector *b = (Vector *)malloc(sizeof(Vector));
    Vector *c = (Vector *)malloc(sizeof(Vector));
    
    a->x = 0.0;
    a->y = 0.0;
    a->z = 10.0;
    
    
    b->x = 200.0;
    b->y = 0.0;
    b->z = 1.0;
    
    
    c->x = 70.0;
    c->y = 120.0;
    c->z = 1.0;
    
    Matrix *m = compute_converter_barycentric(a, b, c);
    
    /*generate some image*/
    unsigned char *image = malloc(IMAGE_WIDTH * IMAGE_HEIGHT * 4);
    unsigned x, y;
    
    Vector point, *point_baryocentric;
    int pixel_power;
    float pixel_z = 0;
    point.z = 1.0;
    for(x = 0; x < IMAGE_WIDTH; x++)
    {
        for(y = 0; y < IMAGE_HEIGHT; y++)
        {
            pixel_power = 0;
            for (int i = 1; i < 5; i++) {
                for (int j = 1; j < 5; j++) {
                    point.x = (float)x + i*0.20;
                    point.y = (float)y + j*0.20;
                    
                    point_baryocentric = affine_transformation(m, &point);
                    point_baryocentric->z = 1.0 - point_baryocentric->x - point_baryocentric->y;
                    
                    if (point_baryocentric->x >= 0.0 &&
                        point_baryocentric->y >= 0.0 &&
                        point_baryocentric->z >= 0.0) {
                        pixel_power++;
                    }
                }
            }
            
            if (pixel_power) {
                point.x = (float)x + 0.5;
                point.y = (float)y + 0.5;
                point_baryocentric = coord_to_bar(m, &point);
                pixel_z = a->z*point_baryocentric->x + b->z*point_baryocentric->y + c->z*point_baryocentric->z;
                add_fragment_to_buffer(f_buffer, x, y, pixel_z, pixel_power);
            }
            
            
        }
    }
    
    for(x = 0; x < IMAGE_WIDTH; x++)
    {
        for(y = 0; y < IMAGE_HEIGHT; y++)
        {
            Fragment_Container frag_cont = f_buffer->fragment_array[x][y];
            int pixel_power = frag_cont.top_fragment.power;
            
            
            image[4 * IMAGE_WIDTH * y + 4 * x + 0] = (unsigned char)((int) (255.0*(float)pixel_power)/16.0);
            image[4 * IMAGE_WIDTH * y + 4 * x + 1] = (unsigned char)((int) (255.0*(float)pixel_power)/16.0);
            image[4 * IMAGE_WIDTH * y + 4 * x + 2] = (unsigned char)((int) (255.0*(float)pixel_power)/16.0);
            image[4 * IMAGE_WIDTH * y + 4 * x + 3] = 255;
        }
    }
    
    /*run an example*/
    encodeOneStep(filename, image, IMAGE_WIDTH, IMAGE_HEIGHT);
    
    free(image);
    
    return 0;
}

