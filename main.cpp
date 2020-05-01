#include "pr.hpp"
#include <iostream>
#include <memory>
#include <functional>
#include <algorithm>

template <int C>
int nextPowerOf(int M)
{
    int a = 1;
    while(a < M)
    {
        a *= C;
    }
    return a;
}

class PJSequence
{
public:
    void clear()
    {
        _samples.clear();
        _random = pr::Xoshiro128StarStar();
    }
    const glm::vec2* samples() const 
    {
        return _samples.data();
    }
    int size() const 
    {
        return _samples.size();
    }

    void extend(int M)
    {
        int N = _samples.size();
        if (N == 0)
        {
            _samples.emplace_back(
                _random.uniformf(),
                _random.uniformf()
            );
            N = 1;
        }
        _samples.resize(nextPowerOf<4>(M));

        // number of cells
        // N = 1, 4, 16, 64, 256...
        while (N < M)
        {
            extendSequence(N, _samples, _random);
            // printf("generated: %d -> %d\n", N, N * 4);
            N = N * 4;
        }
    }
private:
    /*
        i, j : cell index
        xhalf, yhalf: 0 or 1, these indicate the sub-cell
        n : rows, cols count

        (i,j) cell and sub cell
        +------------------+------------------+
        |(xhalf=0, yhalf=0)|(xhalf=1, yhalf=0)|
        +------------------+------------------+
        |(xhalf=0, yhalf=0)|(xhalf=1, yhalf=1)|
        +------------------+------------------+
    */
    static glm::vec2 generateSamplePoint(int i, int j, int xhalf, int yhalf, int n, pr::IRandomNumberGenerator& random)
    {
        return {
            (i + 0.5f * (xhalf + random.uniformf())) / n,
            (j + 0.5f * (yhalf + random.uniformf())) / n,
        };
    }

    // N: generated sample count
    // samples: sample sequence
    static void extendSequence(int N, std::vector<glm::vec2> &samples, pr::IRandomNumberGenerator &random)
    {
        // number of rows, cols
        // n = 1, 2, 4, 8, 16
        int n = sqrt(N);
        for (int s = 0; s < N; ++s)
        {
            glm::vec2 oldpt = samples[s];
            int i = n * oldpt.x;
            int j = n * oldpt.y;
            i = glm::clamp(i, 0, n - 1);
            j = glm::clamp(j, 0, n - 1);

            // (n * oldpt.x - i), (n * oldpt.y - j) : [0, 1) value, normalized coordinates against parent cell
            // 2.0f * (n * oldpt.x - i), 2.0f * (n * oldpt.y - j) : [0, 2) value. these indicates children cell index
            int xhalf = 2.0f * (n * oldpt.x - i);
            int yhalf = 2.0f * (n * oldpt.y - j);
            xhalf = glm::clamp(xhalf, 0, 1);
            yhalf = glm::clamp(yhalf, 0, 1);

            /* choose a diagonal child cell
            +-+-+
            |o| |
            +-+-+
            | |x|
            +-+-+

            o: first cell
            x: diagonal cell
            */
            xhalf = 1 - xhalf;
            yhalf = 1 - yhalf;
            samples[N + s] = generateSamplePoint(i, j, xhalf, yhalf, n, random);

            /* choose a or b
            +-+-+
            |o|a|
            +-+-+
            |b|o|
            +-+-+
            */
            if (random.uniformf() < 0.5f)
            {
                xhalf = 1 - xhalf;
            }
            else
            {
                yhalf = 1 - yhalf;
            }
            samples[2 * N + s] = generateSamplePoint(i, j, xhalf, yhalf, n, random);

            /* choose the last one against the previous p
            +-+-+
            |o|x|
            +-+-+
            |p|o|
            +-+-+
            */
            xhalf = 1 - xhalf;
            yhalf = 1 - yhalf;
            samples[3 * N + s] = generateSamplePoint(i, j, xhalf, yhalf, n, random);
        }
    }
private:
    pr::Xoshiro128StarStar _random;
    std::vector<glm::vec2> _samples;
};

int main() {
    using namespace pr;

    int numberOfSample = 512;

    pr::Xoshiro128StarStar random;
    PJSequence sequence;
    sequence.extend(numberOfSample);
    
    int drawCount = numberOfSample;

    Config config;
    config.ScreenWidth = 1920;
    config.ScreenHeight = 1080;
    config.SwapInterval = 1;
    Initialize(config);

    Camera3D camera;
    camera.origin = { 0, 0, 4 };
    camera.lookat = { 0, 0, 0 };
    camera.zUp = false;

    double e = GetElapsedTime();

    while (pr::NextFrame() == false) {
        if (IsImGuiUsingMouse() == false) {
            UpdateCameraBlenderLike(&camera);
        }

        ClearBackground(0.1f, 0.1f, 0.1f, 1);

        BeginCamera(camera);

        PushGraphicState();

        DrawGrid(GridAxis::XY, 1.0f, 10, { 128, 128, 128 });
        // DrawXYZAxis(1.0f);

        for (int i = 0; i < drawCount; ++i)
        {
            auto x = sequence.samples()[i].x;
            auto y = sequence.samples()[i].y;
            DrawPoint({ x, y, 0.0f }, { 255, 255, 0 }, 4);

            // Projection X
            DrawPoint({ x, 0.0f, 0.0f }, { 255, 0, 0 }, 3);

            // Projection Y
            DrawPoint({ 0.0f, y, 0.0f }, { 0, 255, 0 }, 3);
        }

        PopGraphicState();
        EndCamera();

        BeginImGui();

        ImGui::SetNextWindowSize({ 500, 800 }, ImGuiCond_Once);
        ImGui::Begin("Panel");
        ImGui::Text("fps = %f", GetFrameRate());
        ImGui::SliderInt("draw count", &drawCount, 0, numberOfSample);

        ImGui::End();

        EndImGui();
    }

    pr::CleanUp();
}
