//Modified from http://www.cse.ohio-state.edu/~hwshen/5542/Site/Slides_files/cubeXform.frag
//and http://www.cse.ohio-state.edu/~hwshen/5542/Site/Slides_files/shading_glsl.pdf
//and based on the code presented by graduate teaching associate Soumya Dutta in real time rendering class

//This shader implements the equations for local illumination

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

varying vec4 pcolor; 
varying vec3 v_normal;
varying vec4 pos_in_eye; 

void main() {
	//Prepare vectors
	vec3 lightVector = normalize(vec3(lightPosition - pos_in_eye));
	vec3 eyeVector = normalize(vec3(eyePosition - pos_in_eye));
	vec3 normal = normalize(vec3(v_normal));

	//Ambient
	vec4 ambient = ambient_coef * lightAmbient;
	vec4 ambientBack = ambient_back_coef * lightAmbient;

	//Diffuse
	float ndotl = max(dot(normal,lightVector),0);
	vec4 diffuse = diffuse_coef * lightDiffuse * ndotl;
	float ndot2 = max(dot(-normal,lightVector),0);
	vec4 diffuseBack = diffuse_back_coef * lightDiffuse * ndot2;

	//Specular
	vec3 reflectVector = normalize(reflect(lightVector, normal));
	float rdote = max(dot(reflectVector, eyeVector),0);
	vec4 specular = specular_coef * lightSpecular * pow(rdote, mat_shininess);
	vec3 reflectVectorBack = normalize(reflect(lightVector, -normal));
	float rdoteBack = max(dot(reflectVectorBack, eyeVector),0);
	vec4 specularBack = specular_coef * lightSpecular * pow(rdoteBack, mat_shininess);
	
	//Set the color
	vec4 finalColor = (ambient + diffuse + specular) * pcolor;
	vec4 finalColorBack = (ambientBack + diffuseBack + specularBack) * pcolor;

	if (gl_FrontFacing)
	{
		gl_FragColor = finalColor; 
	}
	else
	{
		gl_FragColor = finalColorBack; 
	}
    
	
	
	//gl_FragColor = vec4(1.0, 1.0, 1.0, 1.0);

 } 