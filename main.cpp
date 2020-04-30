#include "pr.hpp"
#include <iostream>
#include <memory>

pr::Xoshiro128StarStar random;

int N = 0;
glm::vec2 samples[1024];

// i, j : cell index
// xhalf, yhalf: 0 or 1, these indicate the cell
// n : rows, cols
glm::vec2 generateSamplePoint(int i, int j, int xhalf, int yhalf, int n)
{
    return {
        (i + 0.5f * (xhalf + random.uniformf())) / n,
        (j + 0.5f * (yhalf + random.uniformf())) / n,
    };
}

void extendSequence(int N)
{
    // [0]: [0] [1, 2, 3]
    // [0, 1, 2, 3]: [4,5,6], [7,8,9], [10,11,12], [13,14,15]

    // number of rows, cols
    // n = 1, 2, 4, 8, 16
    int n = sqrt(N);
    for (int s = 0; s < N; ++s)
    {
        glm::vec2 oldpt = samples[s];
        int i = n * oldpt.x;
        int j = n * oldpt.y;

        // (n * oldpt.x - i), (n * oldpt.y - j) : [0, 1) value, normalized coordinates against parent cell
        // 2.0f * (n * oldpt.x - i), 2.0f * (n * oldpt.y - j) : [0, 2) value. these indicates children cell index
        int xhalf = 2.0f * (n * oldpt.x - i);
        int yhalf = 2.0f * (n * oldpt.y - j);

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
        samples[N + s] = generateSamplePoint(i, j, xhalf, yhalf, n);

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
        samples[2 * N + s] = generateSamplePoint(i, j, xhalf, yhalf, n);

        /* choose the last one against the previous p
        +-+-+
        |o|x|
        +-+-+
        |p|o|
        +-+-+
        */
        xhalf = 1 - xhalf;
        yhalf = 1 - yhalf;
        samples[3 * N + s] = generateSamplePoint(i, j, xhalf, yhalf, n);
    }
}

void generatePJ(int M)
{
    samples[0] = glm::vec2(random.uniformf(), random.uniformf());
    N = 1;
    // number of cells
    // N = 1, 4, 16, 64, 256...
    while (N < M)
    {
        extendSequence(N);
        printf("generated: %d -> %d\n", N, N * 4);
        N = N * 4;
    }
}

int main() {
    using namespace pr;

    int gen = 512;
    generatePJ(gen);

    int drawCount = gen;

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
        DrawXYZAxis(1.0f);

        for (int i = 0; i < drawCount; ++i)
        {
            DrawCircle({ samples[i].x, samples[i].y, 0.0f }, { 0, 0, 1 }, { 255, 255, 0 }, 0.01f);
        }

        PopGraphicState();
        EndCamera();

        BeginImGui();

        ImGui::SetNextWindowSize({ 500, 800 }, ImGuiCond_Once);
        ImGui::Begin("Panel");
        ImGui::Text("fps = %f", GetFrameRate());
        ImGui::SliderInt("draw count", &drawCount, 0, gen);

        ImGui::End();

        EndImGui();
    }

    pr::CleanUp();
}
