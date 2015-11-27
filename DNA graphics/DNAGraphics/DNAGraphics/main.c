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

int max(int, int);
int min(int, int);

typedef struct matrix3 {
    float a;
    float b;
    float c;
    float d;
    float e;
    float f;
    float g;
    float h;
    float i;
} matrix3;

typedef struct affine_matrix3 {
    float a;
    float b;
    float c;
    float d;
    float tx;
    float ty;
} affine_matrix3;

typedef struct affine_matrix4 {
    float a;
    float b;
    float c;
    float d;
    float e;
    float f;
    float g;
    float h;
    float i;
    float tx;
    float ty;
    float tz;
} affine_matrix4;

typedef struct vector3 {
    float x;
    float y;
    float z;
} vector3;

typedef struct vector4 {
    float x;
    float y;
    float z;
    float w;
} vector4;

typedef struct color {
    int red;
    int green;
    int blue;
} color;

typedef struct vertex {
    vector3 position;
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
    color clr;
} Fragment;

typedef struct fragment_container {
    Fragment top_fragment;
    Fragment *second_frag;
} Fragment_Container;

typedef struct fragment_buffer {
    Fragment_Container fragment_array[IMAGE_WIDTH][IMAGE_HEIGHT];
} fragment_buffer;

//Convenience methods

color * create_color(int red, int green, int blue) {
    color *new_clr = (color *)malloc(sizeof(color));
    new_clr->red = red;
    new_clr->green = green;
    new_clr->blue = blue;
    return(new_clr);
}

color * create_white() {
    color *new_clr = (color *)malloc(sizeof(color));
    new_clr->red = 255;
    new_clr->green = 255;
    new_clr->blue = 255;
    return(new_clr);
}

color * create_red() {
    color *new_clr = (color *)malloc(sizeof(color));
    new_clr->red = 255;
    new_clr->green = 0;
    new_clr->blue = 0;
    return(new_clr);
}

color * create_green() {
    color *new_clr = (color *)malloc(sizeof(color));
    new_clr->red = 0;
    new_clr->green = 255;
    new_clr->blue = 0;
    return(new_clr);
}

color * create_blue() {
    color *new_clr = (color *)malloc(sizeof(color));
    new_clr->red = 0;
    new_clr->green = 0;
    new_clr->blue = 255;
    return(new_clr);
}

vector3 * create_vector3(float x, float y, float z) {
    vector3 *new_pst = (vector3 *)malloc(sizeof(vector3));
    new_pst->x = x;
    new_pst->y = y;
    new_pst->z = z;
    return(new_pst);
}

vector4 * create_vector4(float x, float y, float z, float w) {
    vector4 *new_vec = (vector4 *)malloc(sizeof(vector4));
    new_vec->x = x;
    new_vec->y = y;
    new_vec->z = z;
    new_vec->w = w;
    return(new_vec);
}

triangle * create_triangle(vertex *a, vertex *b, vertex *c) {
    triangle *new_tr = (triangle *)malloc(sizeof(triangle));
    new_tr->vertex_a = *a;
    new_tr->vertex_b = *b;
    new_tr->vertex_c = *c;
    return new_tr;
}

vertex * create_vertex(vector3 *position, color *color) {
    vertex *new_vtx = (vertex *)malloc(sizeof(vertex));
    new_vtx->position = *position;
    new_vtx->clr = *color;
    return new_vtx;
}

matrix3 * create_unit_matrix3() {
    matrix3 *m = (matrix3 *)malloc(sizeof(matrix3));
    m->a = 1.0;
    m->e = 1.0;
    m->i = 1.0;
    return m;
}

affine_matrix3 * create_affine_matrix3(float a, float b, float c, float d, float tx, float ty) {
    affine_matrix3 *m = (affine_matrix3 *)malloc(sizeof(affine_matrix3));
    m->a = a;
    m->b = b;
    m->c = c;
    m->d = d;
    m->tx = tx;
    m->ty = ty;
    return m;
}

affine_matrix3 * create_unit_affine3() {
    //Implement as a macro!
    return create_affine_matrix3(1.0, 0.0, 0.0, 1.0, 0.0, 0.0);
}

affine_matrix3 * create_affine_translation3(int tx, int ty) {
    //Implement as a macro!
    return create_affine_matrix3(1.0, 0.0, 0.0, 1.0, tx, ty);
}

affine_matrix4 * create_affine_matrix4(float a, float b, float c, float d, float e, float f, float g, float h, float i, float tx, float ty, float tz) {
    affine_matrix4 *m = (affine_matrix4 *)malloc(sizeof(affine_matrix4));
    m->a = a;
    m->b = b;
    m->c = c;
    m->d = d;
    m->e = e;
    m->f = f;
    m->g = g;
    m->h = h;
    m->i = i;
    m->tx = tx;
    m->ty = ty;
    m->tz = tz;
    return m;
}

affine_matrix4 * create_unit_affine4() {
    //Implement as a macro!
    return create_affine_matrix4(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
}

affine_matrix4 * create_affine_translation4(int tx, int ty, int tz) {
    //Implement as a macro!
    return create_affine_matrix4(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, tx, ty, tz);
}

vector3 * matric3_mult_vector(matrix3 *m, vector3 *v) {
    vector3 *r = (vector3 *)malloc(sizeof(vector3));
    r->x = m->a*v->x + m->d*v->y + m->g*v->z;
    r->y = m->b*v->x + m->e*v->y + m->h*v->z;
    r->z = m->c*v->x + m->f*v->y + m->i*v->z;
    return r;
}

vector3 * affine_transformation3(affine_matrix3 *m, vector3 *v) {
    vector3 *r = (vector3 *)malloc(sizeof(vector3));
    r->x = m->a*v->x + m->c*v->y + m->tx*v->z;
    r->y = m->b*v->x + m->d*v->y + m->ty*v->z;
    r->z = v->z;
    return r;
}

vector4 * affine_transformation4(affine_matrix4 *m, vector4 *v) {
    vector4 *r = (vector4 *)malloc(sizeof(vector4));
    r->x = m->a*v->x + m->d*v->y + m->g*v->z + m->tx*v->w;
    r->y = m->b*v->x + m->e*v->y + m->h*v->z + m->ty*v->w;
    r->z = m->c*v->x + m->f*v->y + m->i*v->z + m->tz*v->w;
    r->w = v->w;
    return r;
}

vector3 * coord_to_bar(matrix3 *m, vector3 *v) {
    vector3 *r = (vector3 *)malloc(sizeof(vector3));
    r->x = m->a*v->x + m->d*v->y + m->g*v->z;
    r->y = m->b*v->x + m->e*v->y + m->h*v->z;
    r->z = 1 - r->x - r->y;
    return r;
}

matrix3 * matrix3_mult_matrix3(matrix3 *l, matrix3 *r) {
    matrix3 *m = (matrix3 *)malloc(sizeof(matrix3));
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

affine_matrix3 * affine3_mult_affine(affine_matrix3 *l, affine_matrix3 *r) {
    affine_matrix3 *m = (affine_matrix3 *)malloc(sizeof(affine_matrix3));
    return m;
}

void print_vector2(vector3 *v) {
    printf("%f\n%f\n%f\n", v->x, v->y, v->z); //Temp!
}

void print_vector3(vector3 *v) {
    printf("%f\n%f\n%f\n", v->x, v->y, v->z);
}

void print_vector4(vector4 *v) {
    printf("%f\n%f\n%f\n%f\n", v->x, v->y, v->z, v->w);
}

void print_color(color *clr) {
    printf("Red: %d\nGreen: %d\nBlue: %d\n", clr->red, clr->green, clr->blue);
}

void print_matrix3(matrix3 *m) {
    printf("%f %f %f\n%f %f %f\n%f %f %f\n", m->a, m->d, m->g, m->b, m->e, m->h, m->c, m->f, m->i);
}

void print_affine_matrix3(affine_matrix3 *m) {
    printf("%f %f %f\n%f %f %f\n%f %f %f\n", m->a, m->c, m->tx, m->b, m->d, m->ty, 0.0, 0.0, 1.0);
}

void print_affine_matrix4(affine_matrix4 *m) {
    printf("%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n", m->a, m->d, m->g, m->tx, m->b, m->e, m->h, m->ty, m->c, m->f, m->i, m->tz, 0.0, 0.0, 0.0, 1.0);
}

matrix3 * compute_converter_barycentric(vector3 *a, vector3 *b, vector3 *c) {
    matrix3 *m = (matrix3 *)malloc(sizeof(matrix3));
    matrix3 *t = create_unit_matrix3();
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
    return matrix3_mult_matrix3(t, m);
}

color * color_mult(color *color_a, color *color_b, color *color_c, float frac_a, float frac_b, float frac_c) {
    int new_red = (int)(color_a->red*frac_a + color_b->red*frac_b + color_c->red*frac_c);
    int new_green = (int)(color_a->green*frac_a + color_b->green*frac_b + color_c->green*frac_c);
    int new_blue = (int)(color_a->blue*frac_a + color_b->blue*frac_b + color_c->blue*frac_c);
    
    new_red = max(new_red, 0);
    new_green = max(new_green, 0);
    new_blue = max(new_blue, 0);
    new_red = min(new_red, 255);
    new_green = min(new_green, 255);
    new_blue = min(new_blue, 255);
    
    color *new_clr = create_color(new_red, new_green, new_blue);
    
    return new_clr;
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
            return a;
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

fragment_buffer * new_fragment_buffer() {
    fragment_buffer *new_buffer = (fragment_buffer *)malloc(sizeof(fragment_buffer) + sizeof(Fragment_Container)*(IMAGE_WIDTH*IMAGE_HEIGHT - 1));
    return new_buffer;
}

int add_fragment_to_buffer(fragment_buffer *buffer, int x, int y, float z, int power, color *clr) {
    //Extra links don't get dealloced!
    
    Fragment_Container curr_cont;
    curr_cont = buffer->fragment_array[x][y];
    
    //printf("z: %f\n", curr_cont.top_fragment.z);
    if (z < curr_cont.top_fragment.z) {
        if (power == 16) {
            
            //Overrite pixel completely
            curr_cont.top_fragment.z = z;
            curr_cont.top_fragment.power = power;
            curr_cont.top_fragment.clr = *clr;
            //printf("Color red: %d\n", curr_cont.top_fragment.clr.red);
            curr_cont.second_frag = NULL;
            buffer->fragment_array[x][y] = curr_cont;
            return 1;
            //Add macro!
        } else {
            //partially overwrite pixel
            Fragment *old_frag = (Fragment *)malloc(sizeof(Fragment));
            *old_frag = curr_cont.top_fragment;
            curr_cont.top_fragment.z = z;
            curr_cont.top_fragment.power = power;
            curr_cont.top_fragment.clr = *clr;
            buffer->fragment_array[x][y] = curr_cont;
            return 1; //Add macro!
        }
    } else {
        if (curr_cont.top_fragment.power == 16) {
            return 0;
        } else {
            //Try adding fragment to extra link
            
            if (!curr_cont.second_frag || curr_cont.second_frag->z > z) {
                if (!curr_cont.second_frag) {
                    curr_cont.second_frag = (Fragment *)malloc(sizeof(Fragment));
                }
                curr_cont.second_frag->z = z;
                curr_cont.second_frag->power = power;
                curr_cont.second_frag->clr = *clr;
                buffer->fragment_array[x][y] = curr_cont;
                return 1; //Add extra macro!
            } else {
                //Fragment discarded
                return 0;
            }
        }
    }
}

int process_triangle(fragment_buffer *buffer, triangle *triangle) {
    int min_x = (int)float_min(triangle->vertex_a.position.x, triangle->vertex_b.position.x, triangle->vertex_c.position.x);
    int min_y = (int)float_min(triangle->vertex_a.position.y, triangle->vertex_b.position.y, triangle->vertex_c.position.y);
    int max_x = (int)float_max(triangle->vertex_a.position.x, triangle->vertex_b.position.x, triangle->vertex_c.position.x) + 1;
    int max_y = (int)float_max(triangle->vertex_a.position.y, triangle->vertex_b.position.y, triangle->vertex_c.position.y) + 1;

    if (min(min_x, min_y) < 0 || max_x - 1 > IMAGE_HEIGHT || max_y - 1 > IMAGE_HEIGHT) {
        printf("Triangle out of range!\n");
        printf("max_x = %d", max_x);
        return 0;
    }
    
    matrix3 *m = compute_converter_barycentric(&triangle->vertex_a.position, &triangle->vertex_b.position, &triangle->vertex_c.position);
    
    vector3 point, *point_baryocentric;
    int pixel_power;
    float pixel_z = 0;
    color pixel_color;
    point.z = 1.0;
    for(int x = min_x; x < max_x; x++)
    {
        for(int y = min_y; y < max_y; y++)
        {
            pixel_power = 0;
            for (int i = 1; i < 5; i++) {
                for (int j = 1; j < 5; j++) {
                    point.x = (float)x + i*0.20;
                    point.y = (float)y + j*0.20;
                    
                    point_baryocentric = matric3_mult_vector(m, &point);
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
                pixel_z = triangle->vertex_a.position.z*point_baryocentric->x + triangle->vertex_b.position.z*point_baryocentric->y + triangle->vertex_c.position.z*point_baryocentric->z;
                
                pixel_color = *color_mult(&triangle->vertex_a.clr, &triangle->vertex_b.clr, &triangle->vertex_c.clr, point_baryocentric->x, point_baryocentric->y, point_baryocentric->z);
                
                add_fragment_to_buffer(buffer, x, y, pixel_z, pixel_power, &pixel_color);
            }

        }
    }
    return 1; //Processed succesfully
}

int fill_background(fragment_buffer *buffer) {
    for(int x = 0; x < IMAGE_WIDTH; x++)
    {
        for(int y = 0; y < IMAGE_HEIGHT; y++)
        {
            color *black = create_color(0, 0, 0);
            Fragment_Container frag_curr;
            frag_curr.second_frag = NULL;
            frag_curr.top_fragment.z = 10000.0; //Replace by macro!
            frag_curr.top_fragment.clr = *black;
            frag_curr.top_fragment.power = 16;
            buffer->fragment_array[x][y] = frag_curr;
        }
    }
    
    return 1;
}

unsigned char * fill_image_array(fragment_buffer *f_buffer) {
    unsigned char *image = malloc(IMAGE_WIDTH * IMAGE_HEIGHT * 4);
    
    color pixel_color;
    for(int x = 0; x < IMAGE_WIDTH; x++)
    {
        for(int y = 0; y < IMAGE_HEIGHT; y++)
        {
            Fragment_Container frag_cont = f_buffer->fragment_array[x][y];
            if (frag_cont.top_fragment.power == 16 || frag_cont.second_frag == NULL) { //Add macro!
                pixel_color = frag_cont.top_fragment.clr;
                //printf("red: %d\n", pixel_color.red);
            } else {
                int top_power = frag_cont.top_fragment.power;
                int bottom_power = frag_cont.second_frag->power;
                if (top_power + bottom_power <= 16) {
                    pixel_color = *color_mult(&frag_cont.top_fragment.clr, &frag_cont.second_frag->clr, create_white(), ((float)top_power)/16.0, ((float)bottom_power)/16.0, 0.0);
                } else {
                    pixel_color = *color_mult(&frag_cont.top_fragment.clr, &frag_cont.second_frag->clr, create_white(), ((float)top_power)/16.0, ((float)(16 - top_power))/16.0, 0.0);
                }
            }
            
            
            
            image[4 * IMAGE_WIDTH * y + 4 * x + 0] = (unsigned char)((int) ((float)pixel_color.red));
            image[4 * IMAGE_WIDTH * y + 4 * x + 1] = (unsigned char)((int) ((float)pixel_color.green));
            image[4 * IMAGE_WIDTH * y + 4 * x + 2] = (unsigned char)((int) ((float)pixel_color.blue));
            image[4 * IMAGE_WIDTH * y + 4 * x + 3] = 255;
        }
    }
    return image;
}

void encodeOneStep(const char* filename, const unsigned char* image, unsigned width, unsigned height)
{
    /*Encode the image*/
    unsigned int error = lodepng_encode32_file(filename, image, width, height);
    
    /*if there's an error, display it*/
    if (error) printf("error %u: %s\n", error, lodepng_error_text(error));
}

int main(int argc, const char * argv[]) {
    printf("Program terminating!\n");
    fragment_buffer *f_buffer = new_fragment_buffer();

    vector3 *a_pos = create_vector3(0.0, 0.0, 10.0);
    vector3 *b_pos = create_vector3(200.0, 0.0, 10.0);
    vector3 *c_pos = create_vector3(100.0, 200.0, 10.0);
    
    vertex *a = create_vertex(a_pos, create_red());
    vertex *b = create_vertex(b_pos, create_green());
    vertex *c = create_vertex(c_pos, create_blue());

    triangle *tr = create_triangle(a, b, c);

    vertex *d = create_vertex(create_vector3(100.0, 30.0, 1.0), create_white());
    vertex *e = create_vertex(create_vector3(140.0, 60.0, 1.0), create_blue());
    vertex *f = create_vertex(create_vector3(80.0, 130.0, 1.0), create_color(255, 0.0, 1.0));
    
    triangle *tr2 = create_triangle(d, e, f);
    
    fill_background(f_buffer);
    process_triangle(f_buffer, tr2);
    
    process_triangle(f_buffer, tr);
    unsigned char *image = fill_image_array(f_buffer);
    
    //create image
    encodeOneStep(filename, image, IMAGE_WIDTH, IMAGE_HEIGHT);
    
    free(image); //Destroy array
    
    return 0;
}

