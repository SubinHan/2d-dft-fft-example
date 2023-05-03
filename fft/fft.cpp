#include <iostream>

struct Complex
{
    Complex() : real(0.0f), imag(0.0f) {}
    Complex(float real, float imag) : real(real), imag(imag) {}

    float real;
    float imag;

    Complex operator*(const Complex& rhs) const
    {
        return Complex{ real * rhs.real - imag * rhs.imag, real * rhs.imag + imag * rhs.real };
    }

    Complex operator/(const int& rhs) const
    {
        return Complex{ real / rhs, imag / rhs };
    }

    Complex operator+(const Complex& rhs) const
    {
        return Complex{ real + rhs.real, imag + rhs.imag };
    }

    Complex operator-(const Complex& rhs) const
    {
        return Complex{ real - rhs.real, imag - rhs.imag };
    }
};

void dft1d();
void dft2d();
void fft1d();
void fft2d();
void fft(Complex* a, const int size, bool inverse);

int main()
{
    fft2d();
}

void dft1d()
{
    constexpr float PI = 3.141592f;
    constexpr int N = 8;
    const float dt = 0.02f;
    const float fs = 1.0f / dt;
    const float T = 0.16;
    float t[N];

    for(int i = 0; i < N; ++i)
    {
        t[i] = static_cast<float>(i);
    }

    float y[N];
    for(int i = 0; i < N; ++i)
    {
        y[i] = 5.0f + cosf(2.0f * PI * 12.5f * t[i]) + sinf(2.0f * PI * 18.75f * t[i]);
    }

    std::cout << "points:" << std::endl;
    for(int i = 0 ; i < N; ++i)
    {
        printf("%.2f  ", y[i]);
    }
    std::cout << std::endl;
    Complex result[N];
    

    for(int k = 0; k < N; ++k)
    {
	    for(int n = 0; n < N; n++)
	    {
            result[k].real += y[n] * cosf(2.0f * PI * static_cast<float>(k) * static_cast<float>(n) / static_cast<float>(N));
            result[k].imag += y[n] * sinf(2.0f * PI * static_cast<float>(k) * static_cast<float>(n) / static_cast<float>(N));
	    }
    }

    for(int i = 0; i < N; ++i)
    {
        printf("%.2f, %.2f   ", result[i].real, result[i].imag);
    }
    std::cout << std::endl;


    float inversed[N];
    for (int k = 0; k < N; ++k)
    {
        inversed[k] = 0.0f;
        for (int n = 0; n < N; n++)
        {
            inversed[k] += result[n].real * cosf(2.0f * PI * static_cast<float>(k) * static_cast<float>(n) / static_cast<float>(N));
            inversed[k] -= result[n].imag * sinf(2.0f * PI * static_cast<float>(k) * static_cast<float>(n) / static_cast<float>(N));
        }
        inversed[k] /= static_cast<float>(N);
    }

    for (int i = 0; i < N; ++i)
    {
        std::cout << inversed[i] << "  ";
    }
}

void DFT(float* samples, float* complexR, float* complexI, const int size)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            int ij = i * size + j;
            complexR[ij] = 0.0f;
            complexI[ij] = 0.0f;
        }
    }

	constexpr float PI = 3.141592f;
	for(int k = 0; k < size; ++k)
	{
		for (int l = 0; l < size; ++l)
		{
			for (int m = 0; m < size; ++m)
			{
				for(int n = 0; n < size; ++n)
				{
                    const int kl = k * size + l;
                    const int mn = m * size + n;
					complexR[kl] += samples[mn] * cosf(2.0f * PI * 
						(static_cast<float>(k) * static_cast<float>(m) / static_cast<float>(size) + static_cast<float>(l) * static_cast<float>(n) / static_cast<float>(size))) / static_cast<float>(size * size);
					complexI[kl] += samples[mn] * sinf(2.0f * PI * 
						(static_cast<float>(k) * static_cast<float>(m) / static_cast<float>(size) + static_cast<float>(l) * static_cast<float>(n) / static_cast<float>(size))) / static_cast<float>(size * size);
				}
			}
		}
	}
}

void IDFT(float* complexR, float* complexI, float* result, const int size)
{
	constexpr const float PI = 3.141592f;
	for(int k = 0; k < size; ++k)
	{
		for (int l = 0; l < size; ++l)
		{
			for (int m = 0; m < size; ++m)
			{
				for(int n = 0; n < size; ++n)
				{
                    const int kl = k * size + l;
                    const int mn = m * size + n;

                    result[kl] += complexR[mn] * cosf(2.0f * PI *
                        (static_cast<float>(k) * static_cast<float>(m) / static_cast<float>(size) + static_cast<float>(l) * static_cast<float>(n) / static_cast<float>(size)));
                    result[kl] -= complexI[mn] * sinf(2.0f * PI *
                        (static_cast<float>(k) * static_cast<float>(m) / static_cast<float>(size) + static_cast<float>(l) * static_cast<float>(n) / static_cast<float>(size)));
				}
			}
		}
	}
}

void dft2d()
{
    constexpr float PI = 3.141592f;
    constexpr int N = 8;

    float y[N][N];
    for (int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            y[i][j] = 0.0f;
        }
    }

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            y[i][j] += 5.0f + cosf(2.0f * PI * 9.375f * static_cast<float>(j)) + sinf(2.0f * PI * 3.125f * static_cast<float>(j));
            y[j][i] += 5.0f + cosf(2.0f * PI * 6.25f * static_cast<float>(j)) + sinf(2.0f * PI * 6.25f * static_cast<float>(j));
        }
    }

    std::cout << "points:" << std::endl;
    for (int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            printf("%.2f  ", y[i][j]);
        }
        std::cout << std::endl;
    }

    float complexR[N][N];
    float complexI[N][N];


    DFT(&y[0][0], &complexR[0][0], &complexI[0][0], N);

    for (int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            printf("%.2f, %.2f  ", complexR[i][j], complexI[i][j]);
        }
        std::cout << std::endl;
    }

    float inversed[N][N];

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            inversed[i][j] = 0.0f;
        }
    }

    for (int k = 0; k < N; ++k)
    {
        for (int l = 0; l < N; ++l)
        {
            for (int m = 0; m < N; ++m)
            {
                for (int n = 0; n < N; ++n)
                {
                    inversed[k][l] += complexR[m][n] * cosf(2.0f * PI *
                        (static_cast<float>(k) * static_cast<float>(m) / static_cast<float>(N) + static_cast<float>(l) * static_cast<float>(n) / static_cast<float>(N)));
                    inversed[k][l] -= complexI[m][n] * sinf(2.0f * PI *
                        (static_cast<float>(k) * static_cast<float>(m) / static_cast<float>(N) + static_cast<float>(l) * static_cast<float>(n) / static_cast<float>(N)));
                }
            }
        }
    }


    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            printf("%.2f  ", inversed[i][j]);
        }
        std::cout << std::endl;
    }
}

void fft1d()
{
    constexpr float PI = 3.141592f;
    constexpr int N = 8;

    Complex y[N];

    for (int i = 0; i < N; ++i)
    {
		y[i].real = 5.0f + cosf(2.0f * PI * 12.5f * static_cast<float>(i)) + sinf(2.0f * PI * 18.75f * static_cast<float>(i));
    }

    std::cout << "points:" << std::endl;
    for (int i = 0; i < N; ++i)
    {
    	printf("%.2f  ", y[i].real);
    }
    std::cout << std::endl;

    fft(&y[0], N, false);

    for (int i = 0; i < N; ++i)
    {
    	printf("%.2f, %.2f  ", y[i].real, y[i].imag);
    }
    std::cout << std::endl;


    fft(&y[0], N, true);


    for (int i = 0; i < N; ++i)
    {
        printf("%.2f, %.2f  ", y[i].real, y[i].imag);
    }
    std::cout << std::endl;
}

void fft2d()
{
    constexpr float PI = 3.141592f;
    constexpr int N = 8;
    const float dt = 0.02f;
    const float fs = 1.0f / dt;
    const float T = 0.16;
    float t[N][N];

    Complex y[N][N];

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            y[i][j].real += 5.0f + cosf(2.0f * PI * 9.375f * static_cast<float>(j)) + sinf(2.0f * PI * 3.125f * static_cast<float>(j));
            y[j][i].real += 5.0f + cosf(2.0f * PI * 6.25f * static_cast<float>(j)) + sinf(2.0f * PI * 6.25f * static_cast<float>(j));
        }
    }

    std::cout << "points:" << std::endl;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            printf("%.2f  ", y[i][j].real);
        }
        std::cout << std::endl;
    }

    // fft on each row
    for(int i = 0; i < N; ++i)
    {
        fft(&y[i][0], N, false);
    }

    Complex yTranspose[N][N];
    // transpose
    for(int i = 0; i < N; ++i)
    {
	    for(int j = 0; j < N; ++j)
	    {
            yTranspose[i][j] = y[j][i];
	    }
    }

    // fft on each row(it was column before transpose)
    for(int i = 0; i < N; ++i)
    {
        fft(&yTranspose[i][0], N, false);
    }

    // transpose again
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            y[i][j] = yTranspose[j][i];
        }
    }

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            printf("%.2f, %.2f  ", y[i][j].real, y[i][j].imag);
        }
        std::cout << std::endl;
    }


    // ifft on each row
    for (int i = 0; i < N; ++i)
    {
        fft(&y[i][0], N, true);
    }

    // transpose
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            yTranspose[i][j] = y[j][i];
        }
    }

    // ifft on each row(it was column before transpose)
    for (int i = 0; i < N; ++i)
    {
        fft(&yTranspose[i][0], N, true);
    }

    // transpose again
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            y[i][j] = yTranspose[j][i];
        }
    }

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            printf("%.2f, %.2f  ", y[i][j].real, y[i][j].imag);
        }
        std::cout << std::endl;
    }
}

// a의 크기는 2^N 이어야 하고, IDFT의 경우 inverse = true로 넘김
void fft(Complex* a, const int size, bool inverse) {
    constexpr float PI = 3.141592f;
    int n = size;

    // 메모리 재활용을 위한 정렬
    for (int i = 0; i < n; i++) {
        int rev = 0;
        for (int j = 1, target = i; j < n; j <<= 1, target >>= 1) {
            rev = (rev << 1) + (target & 1);
        }
        if (i < rev) {
	        std::swap(a[i], a[rev]);
        }
    }

    for (int len = 2; len <= n; len <<= 1) {
        float x = 2 * PI / len * (inverse ? -1 : 1);
        Complex diff(cos(x), sin(x));
        for (int i = 0; i < n; i += len) {
            Complex e(1.0f, 0.0f);
            for (int j = 0; j < len / 2; j++) {
                int cur = i + j;
                Complex even = a[cur];
                Complex oddE = a[cur + len / 2] * e;
                a[cur] = even + oddE;
                a[cur + len / 2] = even - oddE;
                e = e * diff;
            }
        }
    }
    if (inverse) {
        for (int i = 0; i < n; i++) {
            a[i] = a[i] / n;
        }
    }
}