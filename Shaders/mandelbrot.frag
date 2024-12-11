
vec2 complexAdd(vec2 a, vec2 b) {
    return a + b;
}

vec2 complexMultiply(vec2 a, vec2 b) {
    return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/iResolution.xy;
    //float scale = 2.0/(iTime*1.0);
    float scale = 3.0;
    vec2 c = vec2(uv.x*scale - 2.0, uv.y*scale - 1.5);
    vec2 z = vec2(0.0,0.0);
    z = iMouse.xy/iResolution.xy - vec2(0.5,0.5);
    vec3 col;
    
    int iter = 20;
    int count = 0;
    
    //loop recursive func
    for(int i=0; i<iter; i++)
    {
        if(length(z)>2.0)
        {
            break; //break if outside set
        }
        z = complexMultiply(z, z);
        z = complexAdd(z, c);
        count++;
    }
    if(count==iter){
        fragColor = vec4(0.0,0.0,0.0,1.0);//color black if in
    }
    else{
        //log smoothing
        float log_zn = log(pow(z.x,2.0) + pow(z.y,2.0)) / 2.0;
        float nu = log(log_zn / log(2.0)) / log(2.0);
        float iteration = float(count) + 1.0 - nu;
        
        
        col = mix(vec3(0.0,0.0,0.0), vec3(1.0,1.0,0.0), iteration/float(iter));
        fragColor = vec4(col, 1.0);
    }
    
}
