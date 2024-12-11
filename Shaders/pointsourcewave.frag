void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/iResolution.xy;
    vec2 mousepos = iMouse.xy/iResolution.xy;
    
    vec2 p1 = vec2(0.2, 0.3);
    vec2 p2 = vec2(mousepos.x, mousepos.y);
        
    float d1 = sqrt(pow(p1.x - uv.x, 2.0) + pow(p1.y - uv.y, 2.0));
    float d2 = sqrt(pow(p2.x - uv.x, 2.0) + pow(p2.y - uv.y, 2.0));
    
    float w1 = 200.0;
    float w2 = 200.0;
    
    float f1 = 40.0;
    float f2 = 40.0;
    
    float wave1 = sin(w1 *d1 - f1*iTime);
    float wave2 = sin(w2 *d2 - f2*iTime);
    
    
    // Output to screen
    //fragColor = vec4(wave1, wave2, wave1+wave2,1.0);
    fragColor = vec4(wave1+wave2, wave1+wave2, wave1+wave2, 1.0);
}