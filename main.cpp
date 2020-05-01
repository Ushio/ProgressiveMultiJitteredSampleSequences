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

class RandomSequence
{
public:
    // this method should be called before first extend. 
    void setSeed(uint32_t s)
    {
        _seed = s;
    }
    void clear() {
        _samples.clear();
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
        int i = _samples.size();
        if (i == 0)
        {
            _random = decltype(_random)(_seed);
        }
        _samples.resize(M);
        for (; i < M; ++i)
        {
            _samples[i] = { _random.uniformf(), _random.uniformf() };
        }
    }
private:
    uint32_t _seed = 1;
    pr::Xoshiro128StarStar _random;
    std::vector<glm::vec2> _samples;
};

class PJSequence
{
public:
    // PJ is not required integer calculation, but we use this for consistency with PMJ
    enum {
        RANDOM_MAX = 0x7FFFFF,
        RANDOM_LENGTH,
    };

    // this method should be called before first extend. 
    void setSeed(uint32_t s)
    {
        _seed = s;
    }
    void clear()
    {
        _samples.clear();
    }
    const glm::ivec2* samples() const
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
            _random = decltype(_random)(_seed);
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
    glm::vec2 to01(glm::ivec2 s) const
    {
        return glm::vec2(s) / glm::vec2(RANDOM_LENGTH);
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
    static glm::ivec2 generateSamplePoint(int i, int j, int xhalf, int yhalf, int n, pr::IRandomNumberGenerator& random)
    {
        int squareLength = (RANDOM_LENGTH / n);
        int halfSquareLength = squareLength / 2;
        int x = i * squareLength + xhalf * halfSquareLength + (random.uniformi() % halfSquareLength);
        int y = j * squareLength + yhalf * halfSquareLength + (random.uniformi() % halfSquareLength);
        return {
            x,
            y,
        };
    }

    // N: generated sample count
    // samples: sample sequence
    static void extendSequence(int N, std::vector<glm::ivec2> &samples, pr::IRandomNumberGenerator &random)
    {
        // number of rows, cols
        // n = 1, 2, 4, 8, 16
        int n = sqrt(N);
        for (int s = 0; s < N; ++s)
        {
            glm::ivec2 oldpt = samples[s];
            int squareLength = (RANDOM_LENGTH / n);
            int i = oldpt.x / squareLength;
            int j = oldpt.y / squareLength;
            int i_mod = oldpt.x % squareLength;
            int j_mod = oldpt.y % squareLength;

            // local sub-square index
            int xhalf = i_mod < (squareLength / 2) ? 0 : 1;
            int yhalf = j_mod < (squareLength / 2) ? 0 : 1;

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

            // b -> a or a -> b
            xhalf = 1 - xhalf;
            yhalf = 1 - yhalf;
            samples[3 * N + s] = generateSamplePoint(i, j, xhalf, yhalf, n, random);
        }
    }
private:
    uint32_t _seed = 1;
    pr::Xoshiro128StarStar _random;
    std::vector<glm::ivec2> _samples;
};

class PMJSequence
{
public:
    // Highly reccommend to use integer for sample. 
    // We have to avoid numerical error because the stratum check must be strict.
    enum {
        RANDOM_MAX = 0x7FFFFF,
        RANDOM_LENGTH,
    };

    // this method should be called before first extend. 
    void setSeed(uint32_t s)
    {
        _seed = s;
    }
    void clear()
    {
        _samples.clear();
    }
    const glm::ivec2* samples() const
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
            _random = decltype(_random)(_seed);
            _samples.emplace_back(
                _random.uniformi() % RANDOM_LENGTH,
                _random.uniformi() % RANDOM_LENGTH
            );
            N = 1;
        }
        _samples.resize(nextPowerOf<4>(M));

        // number of cells
        // N = 1, 4, 16, 64, 256...
        while (N < M)
        {
            std::vector<bool> xstratum, ystratum;
            buildOccupied(N, _samples, xstratum, ystratum);
            extendSequenceDiagonal(N, _samples, _random, xstratum, ystratum);
            // check whether stratums are filled.
            // PR_ASSERT(std::all_of(xstratum.begin(), xstratum.end(), [](bool b) { return b; }), "");
            // PR_ASSERT(std::all_of(ystratum.begin(), ystratum.end(), [](bool b) { return b; }), "");

            buildOccupied(N * 2, _samples, xstratum, ystratum);
            extendSequenceNonDiagonal(N * 2, _samples, _random, xstratum, ystratum);
            // check whether stratums are filled.
            // PR_ASSERT(std::all_of(xstratum.begin(), xstratum.end(), [](bool b) { return b; }), "");
            // PR_ASSERT(std::all_of(ystratum.begin(), ystratum.end(), [](bool b) { return b; }), "");

            // printf("generated: %d -> %d\n", N, N * 4);
            N = N * 4;
        }
    }
    glm::vec2 to01(glm::ivec2 s) const
    {
        return glm::vec2(s) / glm::vec2(RANDOM_LENGTH);
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
    static glm::ivec2 generateSamplePoint(int N, int i, int j, int xhalf, int yhalf, int n, pr::IRandomNumberGenerator& random, std::vector<bool>& xstratum, std::vector<bool>& ystratum)
    {
        /*
         This is the stratum count.
         We'll generate [N, Nx2) samples on the current step.
         So N stratums are already generated and filled so stratums should be Nx2 and these will be filled on the current step.
        */
        int Nx2 = N * 2;

        int squareLength = (RANDOM_LENGTH / n);
        int halfSquareLength = squareLength / 2;
        int stratumLength = (RANDOM_LENGTH / Nx2);
        int x;
        for (;;)
        {
            x = i * squareLength + xhalf * halfSquareLength + (random.uniformi() % halfSquareLength);
            int xstratum_index = x / stratumLength;
            if (xstratum[xstratum_index] == false)
            {
                xstratum[xstratum_index] = true;
                break;
            }
        }

        int y;
        for (;;)
        {
            y = j * squareLength + yhalf * halfSquareLength + (random.uniformi() % halfSquareLength);
            int ystratum_index = y / stratumLength;
            if (ystratum[ystratum_index] == false)
            {
                ystratum[ystratum_index] = true;
                break;
            }
        }

        return {
            x,
            y
        };
    }

    static void buildOccupied(int N, std::vector<glm::ivec2>& samples, std::vector<bool>& xstratum, std::vector<bool>& ystratum)
    {
        /*
         This is the stratum count.
         We'll generate [N, Nx2) samples on the current step.
         So N stratums are already generated and filled so stratums should be Nx2 and these will be filled on the current step.
        */
        int Nx2 = N * 2;

        xstratum.clear();
        ystratum.clear();
        xstratum.resize(Nx2);
        ystratum.resize(Nx2);
        std::fill(xstratum.begin(), xstratum.end(), false);
        std::fill(ystratum.begin(), ystratum.end(), false);

        for (int i = 0; i < N; ++i)
        {
            int xstratum_index = samples[i].x / (RANDOM_LENGTH / Nx2);
            int ystratum_index = samples[i].y / (RANDOM_LENGTH / Nx2);
            xstratum[xstratum_index] = true;
            ystratum[ystratum_index] = true;
        }
    }

    // Generate [N, Nx2) sequence.
    // samples: sample sequence
    static void extendSequenceDiagonal(int N, std::vector<glm::ivec2>& samples, pr::IRandomNumberGenerator& random, std::vector<bool>& xstratum, std::vector<bool>& ystratum)
    {
        // number of rows, cols
        // n = 1, 2, 4, 8, 16
        int n = sqrt(N);
        for (int s = 0; s < N; ++s)
        {
            glm::ivec2 oldpt = samples[s];
            int squareLength = (RANDOM_LENGTH / n);
            int i = oldpt.x / squareLength;
            int j = oldpt.y / squareLength;
            int i_mod = oldpt.x % squareLength;
            int j_mod = oldpt.y % squareLength;

            // local sub-square index
            int xhalf = i_mod < (squareLength / 2) ? 0 : 1;
            int yhalf = j_mod < (squareLength / 2) ? 0 : 1;

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
            samples[N + s] = generateSamplePoint(N, i, j, xhalf, yhalf, n, random, xstratum, ystratum);
        }
    }
    // Generate [Nx2, Nx3) sequence.
    static void extendSequenceNonDiagonal(int Nx2, std::vector<glm::ivec2>& samples, pr::IRandomNumberGenerator& random, std::vector<bool>& xstratum, std::vector<bool>& ystratum)
    {
        int N = Nx2 / 2;
        // number of rows, cols
        // n = 1, 2, 4, 8, 16
        int n = sqrt(N);
        for (int s = 0; s < N; ++s)
        {
            glm::ivec2 oldpt = samples[s];
            int squareLength = (RANDOM_LENGTH / n);
            int i = oldpt.x / squareLength;
            int j = oldpt.y / squareLength;
            int i_mod = oldpt.x % squareLength;
            int j_mod = oldpt.y % squareLength;

            // local sub-square index
            int xhalf = i_mod < (squareLength / 2) ? 0 : 1;
            int yhalf = j_mod < (squareLength / 2) ? 0 : 1;

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

            samples[N * 2 + s] = generateSamplePoint(Nx2, i, j, xhalf, yhalf, n, random, xstratum, ystratum);
 
            // b -> a or a -> b
            xhalf = 1 - xhalf;
            yhalf = 1 - yhalf;
            samples[N * 3 + s] = generateSamplePoint(Nx2, i, j, xhalf, yhalf, n, random, xstratum, ystratum);
        }
    }
private:
    uint32_t _seed = 1;
    pr::Xoshiro128StarStar _random;
    std::vector<glm::ivec2> _samples;
};

enum {
    SAMPLES_RANDOM = 0,
    SAMPLES_PJ,
    SAMPLES_PMJ,
    SAMPLES_TYPE_COUNT
};
const char* Samples[] = {
    "SAMPLES_RANDOM",
    "SAMPLES_PJ",
    "SAMPLES_PMJ",
};

int main() {
    using namespace pr;

    const int numberOfSample = 8096;

    int seed = 1;
    bool autoIncrementSeed = false;
    int drawCount = 128;
    int sampleMode = SAMPLES_PJ;
    int pixelsize = 3;

    RandomSequence random;
    random.setSeed(seed);
    random.extend(numberOfSample);

    PJSequence pj;
    pj.setSeed(seed);
    pj.extend(numberOfSample);

    PMJSequence pmj;
    pmj.setSeed(seed);
    pmj.extend(numberOfSample);

    Config config;
    config.ScreenWidth = 1920;
    config.ScreenHeight = 1080;
    config.SwapInterval = 1;
    Initialize(config);

    Camera3D camera;
    camera.origin = { 0.5f, 0.5f, 2 };
    camera.lookat = { 0.5f, 0.5f, 0 };
    camera.zUp = false;

    double e = GetElapsedTime();

    while (pr::NextFrame() == false) {

        if (IsKeyDown(KEY_DOWN))
        {
            sampleMode = std::min(sampleMode + 1, (int)SAMPLES_TYPE_COUNT - 1);
        }
        if (IsKeyDown(KEY_UP))
        {
            sampleMode = std::max(sampleMode - 1, 0);
        }

        bool seed_update = false;

        if (IsKeyDown(KEY_RIGHT))
        {
            seed += 1;
            seed_update = true;
        }
        if (IsKeyDown(KEY_LEFT))
        {
            seed -= 1;
            seed_update = true;
        }

        if (autoIncrementSeed)
        {
            seed += 1;
            seed_update = true;
        }

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
            float x = 0;
            float y = 0;
            switch (sampleMode) {
            case SAMPLES_RANDOM:
                x = random.samples()[i].x;
                y = random.samples()[i].y;
                break;
            case SAMPLES_PJ: {
                glm::vec2 p = pj.to01(pj.samples()[i]);
                x = p.x;
                y = p.y;
                break;
            }
            case SAMPLES_PMJ: {
                glm::vec2 p = pmj.to01(pmj.samples()[i]);
                x = p.x;
                y = p.y;
                break;
            }
            }

            DrawPoint({ x, y, 0.0f }, { 255, 255, 0 }, pixelsize);

            // Projection X
            DrawPoint({ x, 0.0f, 0.0f }, { 255, 0, 0 }, pixelsize);

            // Projection Y
            DrawPoint({ 0.0f, y, 0.0f }, { 0, 255, 0 }, pixelsize);
        }

        //if(sampleMode != SAMPLES_RANDOM)
        //{
        //    //nextPowerOf<2>(sampleMode)
        //}
        // DrawGrid(GridAxis::XY, 1.0f, 10, { 128, 128, 128 });

        PopGraphicState();
        EndCamera();

        BeginImGui();

        ImGui::SetNextWindowSize({ 500, 800 }, ImGuiCond_Once);
        ImGui::Begin("Panel");
        ImGui::Text("fps = %f", GetFrameRate());

        if (ImGui::SliderInt("seed", &seed, 0, 1024) || seed_update)
        {
            random.setSeed(seed);
            random.clear();
            random.extend(numberOfSample);

            pj.setSeed(seed);
            pj.clear();
            pj.extend(numberOfSample);

            pmj.setSeed(seed);
            pmj.clear();
            pmj.extend(numberOfSample);
        }
        ImGui::Checkbox("auto increment seed", &autoIncrementSeed);
        ImGui::SliderInt("pixel size", &pixelsize, 0, 5);
        ImGui::SliderInt("draw count", &drawCount, 0, numberOfSample);
        ImGui::Combo("Sample Mode", &sampleMode, Samples, IM_ARRAYSIZE(Samples));
        
        ImGui::End();

        EndImGui();
    }

    pr::CleanUp();
}
