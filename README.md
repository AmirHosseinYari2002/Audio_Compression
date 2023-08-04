<h1 align='center'> Audio Compression </h1>

<h2 align='center'> Sharif University of Technology </h2>

<h3 align='center'> Electrical Engineering Department </h3>

MDCT stands for "Modified Discrete Cosine Transform." The MDCT is a modified version of the standard DCT, and it is commonly used in audio compression algorithms like MP3 and AAC. It improves the efficiency of signal representation by splitting the input signal into overlapping blocks before applying the DCT. The overlapping blocks help reduce artifacts caused by the boundary discontinuities that can occur when using a standard DCT on consecutive blocks.

MDCT is a crucial part of many modern audio codecs, as it enables high levels of compression while maintaining acceptable audio quality. By transforming audio data into the frequency domain and quantizing the resulting coefficients, these codecs can significantly reduce the file size of audio recordings without perceptible loss of quality to the listener.

In this project, we have implemented a basic design for audio codec based on MDCT and we have checked the quality of coded sound.
