#define SIGNATURE "Sung-Hun Bae 2009.02"

#include <vector>
#include <cstring>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdlib.h>
#include "gts.h"

#define MAXLENGTH 256

using namespace std;

static void build_list (gpointer data, GSList ** list)
{
  /* always use O(1) g_slist_prepend instead of O(n) g_slist_append */
  *list = g_slist_prepend (*list, data);
}
static void build_list1 (gpointer data, GList ** list)
{
  /* always use O(1) g_list_prepend instead of O(n) g_list_append */
  *list = g_list_prepend (*list, data);
}
static void build_list2 (GtsVertex * v, GList ** list)
{
  if (gts_vertex_is_boundary (v, NULL))
    *list = g_list_prepend (*list, v);
}

static void edge_cleanup (GtsSurface * surface)
{
  GSList * edges = NULL;
  GSList * i;

  g_return_if_fail (surface != NULL);

  /* build list of edges */
  gts_surface_foreach_edge (surface, (GtsFunc) build_list, &edges);

  /* remove degenerate and duplicate edges.
     Note: we could use gts_edges_merge() to remove the duplicates and then
     remove the degenerate edges but it is more efficient to do everything
     at once (and it's more pedagogical too ...) */

  /* We want to control manually the destruction of edges */
  gts_allow_floating_edges = TRUE;

  i = edges;
  while (i) {
    GtsEdge * e = (GtsEdge *) i->data;
    GtsEdge * duplicate;
    if (GTS_SEGMENT (e)->v1 == GTS_SEGMENT (e)->v2) /* edge is degenerate */
      /* destroy e */
      gts_object_destroy (GTS_OBJECT (e));
    else if ((duplicate = gts_edge_is_duplicate (e))) {
      /* replace e with its duplicate */
      gts_edge_replace (e, duplicate);
      /* destroy e */
      gts_object_destroy (GTS_OBJECT (e));
    }
    i = i->next;
  }

  /* don't forget to reset to default */
  gts_allow_floating_edges = FALSE;

  /* free list of edges */
  g_slist_free (edges);
}

static void triangle_cleanup (GtsSurface * s)
{
  GSList * triangles = NULL;
  GSList * i;

  g_return_if_fail (s != NULL);

  /* build list of triangles */
  gts_surface_foreach_face (s, (GtsFunc) build_list, &triangles);

  /* remove duplicate triangles */
  i = triangles;
  while (i) {
    GtsTriangle * t = (GtsTriangle *) i->data;
    if (gts_triangle_is_duplicate (t))
      /* destroy t, its edges (if not used by any other triangle)
         and its corners (if not used by any other edge) */
      gts_object_destroy (GTS_OBJECT (t));
    i = i->next;
  }

  /* free list of triangles */
  g_slist_free (triangles);
}

int main (int argc, char *argv[])
{
  char line[MAXLENGTH],vert[MAXLENGTH],face[MAXLENGTH];
  char gtsfile [MAXLENGTH];
  FILE *msms_vert, *msms_face, *gts;
  unsigned int faces = 0, edges = 0;

  if (argc < 2 || argc > 3) {
    printf("\t\t\t\t\t\t\t" SIGNATURE "\n");
    printf("\tconvert MSMS surface to GTS surface\n\n");
    printf("\tUsage: msms2gts <msms_file_prefix> [<# of faces>]\n");
    printf("\tif # of faces is given, MSMS surfaces are coarsened ");
    printf("to that level\n");
    printf("\n");
    exit(1);
    }

  sprintf(vert,"%s.vert",argv[1]);
  sprintf(face,"%s.face",argv[1]);
  if (argc == 3) 
    sscanf (argv[2],"%u",&faces);
  
  int l,vi=0,va,vb,vc,vnum;
  double vx,vy,vz;
  gsl_matrix *v = NULL;

  GtsSurface *S = gts_surface_new (
		gts_surface_class(),
		gts_face_class(),
		gts_edge_class(),
		gts_vertex_class());

  GtsVertex *v1,*v2,*v3;
  GtsEdge *a,*b,*c;
  GtsFace *dS;

  if ((msms_vert = fopen(vert,"rt")) == NULL) {
    printf("error: cannot open MSMS .vert file\n");
    exit(1);
    }

  if ((msms_face = fopen(face,"rt")) == NULL) {
    printf("error: cannot open MSMS .face file\n");
    exit(1);
    }

  /* read MSMS .vert file */
  l = 0;
  while(fgets(line, sizeof(line), msms_vert)!=NULL) {
    l++;
    if (l == 3) {
      sscanf(line,"%d",&vnum);
      v = gsl_matrix_alloc(vnum,3);
      vi = 0;
      }
    if (l > 3) {
      sscanf(line,"%lf %lf %lf",&vx,&vy,&vz);
      gsl_matrix_set (v,vi,0,vx);
      gsl_matrix_set (v,vi,1,vy);
      gsl_matrix_set (v,vi,2,vz);
      vi++;
      }
    }
  fclose(msms_vert);
	
  /* read MSMS .face file */
  l = 0;
  while(fgets(line, sizeof(line), msms_face)!=NULL) {
    l++;
    if (l > 3) {
      sscanf(line,"%d %d %d",&va,&vb,&vc);
      v1 = gts_vertex_new (S->vertex_class, 
	      gsl_matrix_get(v,va-1,0),
	      gsl_matrix_get(v,va-1,1),
	      gsl_matrix_get(v,va-1,2));
      v2 = gts_vertex_new (S->vertex_class, 
	      gsl_matrix_get(v,vb-1,0),
	      gsl_matrix_get(v,vb-1,1),
	      gsl_matrix_get(v,vb-1,2));
      v3 = gts_vertex_new (S->vertex_class, 
	      gsl_matrix_get(v,vc-1,0),
	      gsl_matrix_get(v,vc-1,1),
	      gsl_matrix_get(v,vc-1,2));
      a  = gts_edge_new (S->edge_class,v1,v2);
      b  = gts_edge_new (S->edge_class,v2,v3);
      c  = gts_edge_new (S->edge_class,v3,v1);
      dS = gts_face_new (S->face_class,a,b,c);
      gts_surface_add_face (S,dS);
      }
    }
  fclose(msms_face);
  gsl_matrix_free (v);

  /* CLEANUP */
  GList * vertices = NULL;
  gdouble threshold = 0.;
  gboolean boundary = FALSE;
  gboolean (* check) (GtsVertex *, GtsVertex *) = NULL;

  /* merge vertices which are close enough */
  /* build list of vertices */
  gts_surface_foreach_vertex (S, boundary ? 
    (GtsFunc) build_list2 : (GtsFunc) build_list1, &vertices);
  /* merge vertices: we MUST update the variable vertices because this 
    function modifies the list (i.e. removes the merged vertices). */
  vertices = gts_vertices_merge (vertices, threshold, check);

  /* free the list */
  g_list_free (vertices);

  /* eliminate degenerate and duplicate edges */
  edge_cleanup (S);

  /* eliminate duplicate triangles */
  triangle_cleanup (S);


  /* COARSENING */
  if (faces > 0) {

    GtsVolumeOptimizedParams params = { 0.5, 0.5, 0. };
    gdouble fold = 3.14159265359/180.;

    edges = faces*3/2;

    gts_surface_coarsen (S, 
	(GtsKeyFunc) gts_volume_optimized_cost, &params,	/* cost */
	(GtsCoarsenFunc) gts_volume_optimized_vertex, &params,  /* coarsen */
	(GtsStopFunc) gts_coarsen_stop_number, &edges,		/* stop */
	fold);
    sprintf(gtsfile,"%s.%04d.gts",argv[1],faces);
    }
  else 
    sprintf(gtsfile,"%s.gts",argv[1]);

  if ((gts = fopen(gtsfile,"wt")) == NULL) {
    printf("error: cannot open GTS .gts file\n");
    exit(1);
    }

  gts_surface_write (S, gts);

  return 0;
}
