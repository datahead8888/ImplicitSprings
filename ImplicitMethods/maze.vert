//Modified from http://www.cse.ohio-state.edu/~hwshen/5542/Site/Slides_files/cubeXform.vert
//and http://www.cse.ohio-state.edu/~hwshen/5542/Site/Slides_files/shading_glsl.pdf
//and based on the code presented by Soumya Dutta in class

//This vertex shader applies the local2clip transformation to each point to set the vertex position, applies the local2eye transformation to each vertex so that the fragment shader can use it,
//transforms normals, and passes info to the fragment shader

attribute vec4 position; 
attribute vec4 color1; 
attribute vec4 normal; 


uniform vec4 ambient_coef;
uniform vec4 diffuse_coef;
uniform vec4 ambient_back_coef;
uniform vec4 diffuse_back_coef;
uniform vec4 specular_coef;
uniform float mat_shininess;
uniform vec4 lightAmbient;
uniform vec4 lightDiffuse;
uniform vec4 lightSpecular;
uniform vec4 lightPosition;
uniform vec4 eyePosition;

uniform mat4 local2clip; 
uniform mat4 local2eye;
uniform mat4 normalMatrix;

varying vec4 pcolor;
varying vec3 v_normal;
varying vec4 pos_in_eye; 

void main(){
	  //Transform the normal using the transpose of the inverse, and transform the vertex into eye space
	  v_normal = normalize(vec3(normalMatrix * normal));
	  pos_in_eye = local2eye * position;
     
      pcolor = color1; 
      
      gl_Position = local2clip * position; 

}