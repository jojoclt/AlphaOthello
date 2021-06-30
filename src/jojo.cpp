#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <queue>
#include <vector>
#include<cassert>
#include<cstring>
#include<cstdint>

#define inf 1e9

const int DEPTH = 5;
const int SIZE = 8;
using ARR = std::array<std::array<int, SIZE>, SIZE>;
using std::vector;

// uint64_t zobristTable[SIZE][SIZE][2];
uint64_t hash = 0;

enum Algo { purerandom, statevalue, minimax, alphabeta, mtdf };
Algo algo = alphabeta;

struct Point {
    int x, y;
    Point() : Point(0, 0) {}
    Point(int x, int y) : x(x), y(y) {}
    bool operator==(const Point& rhs) const { return x == rhs.x && y == rhs.y; }
    bool operator!=(const Point& rhs) const { return !operator==(rhs); }
    Point operator+(const Point& rhs) const
    {
        return Point(x + rhs.x, y + rhs.y);
    }
    Point operator-(const Point& rhs) const
    {
        return Point(x - rhs.x, y - rhs.y);
    }

};

struct Node {
    Point p;
    double val;
    bool operator<(Node t) const { return t.val > val; }
    bool operator>(Node t) const { return t.val < val; }
    Node() : p(0,0), val(0) {}
    Node(Point p, double val) : p(p), val(val) {}
};

int player;
ARR board;
vector<Point> next_valid_spots;
Point bestMove;

enum SPOT_STATE { EMPTY = 0, BLACK = 1, WHITE = 2 };
const std::array<Point, 8> directions{{
        Point(-1, -1), Point(-1, 0), Point(-1, 1),
        Point(0, -1), /*{0, 0}, */Point(0, 1),
        Point(1, -1), Point(1, 0), Point(1, 1)
    }};
const std::array<Point, 4> corners = {{ 
    Point(0, 0), Point(0, 7), Point(7, 0), Point(7, 7) 
}};
int getNextPlayer(int player) { return 3 - player; }
bool is_spot_on_board(Point p){ return 0 <= p.x && p.x < SIZE && 0 <= p.y && p.y < SIZE; }
int get_disc(ARR _board, Point p) { return _board[p.x][p.y]; }
void set_disc(ARR& _board, Point p, int disc) { _board[p.x][p.y] = disc; }
bool is_disc_at(ARR _board, Point p, int disc)
{
    if (!is_spot_on_board(p)) return false;
    if (get_disc(_board, p) != disc) return false;
    return true;
}
bool is_spot_valid(ARR _board, Point center, int curPlayer) {
    if (get_disc(_board, center) != EMPTY)
        return false;
    for (Point dir: directions) {
        // Move along the direction while testing.
        Point p = center + dir;
        if (!is_disc_at(_board, p, getNextPlayer(curPlayer)))
            continue;
        p = p + dir;
        while (is_spot_on_board(p) && get_disc(_board, p) != EMPTY) {
            if (is_disc_at(_board, p, curPlayer))
                return true;
            p = p + dir;
        }
    }
    return false;
}
vector<Point> get_valid_spots(ARR _board, int curPlayer) {
    std::vector<Point> v;
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            Point p = Point(i, j);
            if (_board[i][j] != EMPTY)
                continue;
            if (is_spot_valid(_board, p, curPlayer))
                v.push_back(p);
        }
    }
    return v;
}
void flip_discs(ARR& _board, Point center, int curPlayer)
{
    for (Point dir : directions) {
        // Move along the direction while testing.
        Point p = center + dir;
        if (!is_disc_at(_board, p, getNextPlayer(curPlayer))) continue;
        std::vector<Point> discs({p});
        p = p + dir;
        while (is_spot_on_board(p) && get_disc(_board, p) != EMPTY) {
            if (is_disc_at(_board, p, curPlayer)) {
                for (Point s : discs) {
                    set_disc(_board, s, curPlayer);
                }
                break;
            }
            discs.push_back(p);
            p = p + dir;
        }
    }
}
void put_disc(ARR& _board, Point p, int curPlayer)
{
    set_disc(_board, p, curPlayer);
    flip_discs(_board, p, curPlayer);
}
bool visit[SIZE][SIZE] = {};
int HelperCountVert(ARR _board, Point pos, int dir , int ply) {
    // DIR 1 7 UP DOWN
    int c = 0;
    Point p = pos;
    if (visit[p.x][p.y]) return 0;
    while (_board[p.x][p.y] == ply) {
        c++;
        visit[p.x][p.y] = true;
        p = p + directions[dir];
    }
    return c;
}
int countStableDisc(ARR _board, Point p, int ply) {
    int c = 0;
    int arr[SIZE], i;
    // LOOK DOWN
    if (p == corners[0]) {
        for (i = 0; i < SIZE; i++) {
            arr[i] = HelperCountVert(_board, Point(0,i), 7, ply);
        }
        c += arr[0];
        for (i = 0; i < SIZE - 1; i++) {
            if (arr[i] && arr[i+1] && arr[i] >= arr[i+1]) c += arr[i+1];
        }
    }
    else if (p == corners[1]) {
        for (i = SIZE - 1; i >= 0; i--) {
            arr[i] = HelperCountVert(_board, Point(0,i), 7, ply);
        }
        c += arr[SIZE - 1];
        for (i = SIZE - 1; i >= 0; i--) {
            if (arr[i] && arr[i+1] && arr[i] >= arr[i+1]) c += arr[i+1];
        }
    }

    else if (p == corners[2]) {
         for (i = 0; i < SIZE; i++) {
            arr[i] = HelperCountVert(_board, Point(7,i), 1, ply);
        }
        c += arr[0];
        for (i = 0; i < SIZE - 1; i++) {
            if (arr[i] && arr[i+1] && arr[i] >= arr[i+1]) c += arr[i+1];
        }
    }

    else if (p == corners[3]) {
         for (i = SIZE - 1; i >= 0; i--) {
            arr[i] = HelperCountVert(_board, Point(7,i), 7, ply);
        }
        c += arr[SIZE - 1];
        for (i = SIZE - 1; i >= 0; i--) {
            if (arr[i] && arr[i+1] && arr[i] >= arr[i+1]) c += arr[i+1];
        }
    }
    return c;
}
double Heuristic(ARR _board, int curPlayer)
{
    int count[3] = {};
    double V = 0, D = 0, C = 0, CS = 0, MC = 0, SDC = 0;
    ARR w;
    w[0] = {20, -3, 11, 8, 8, 11, -3, 20};
    w[1] = {-3, -7, -3, -1, -1, -3, -7, -3};
    w[2] = {11, -3, 2, 2, 2, 2, -3, 11};
    w[3] = {8, 1, 2, -3, -3, 2, 1, 8};
    w[4] = {8, 1, 2, -3, -3, 2, 1, 8};
    w[5] = {11, -3, 2, 2, 2, 2, -3, 11};
    w[6] = {-3, -7, -3, -1, -1, -3, -7, -3};
    w[7] = {20, -3, 11, 8, 8, 11, -3, 20};

    // Position Values
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            if (_board[i][j] == curPlayer)
                V += w[i][j];
            else if (_board[i][j] == getNextPlayer(curPlayer))
                V -= w[i][j];
            count[_board[i][j]]++;
        }
    }
    // DISC COUNT
    if (count[curPlayer] > count[getNextPlayer(curPlayer)])
        D = (100.0 * count[curPlayer]) /
            (count[curPlayer] + count[getNextPlayer(curPlayer)]);
    else if (count[curPlayer] < count[getNextPlayer(curPlayer)])
        D = -(100.0 * count[curPlayer]) /
            (count[curPlayer] + count[getNextPlayer(curPlayer)]);
    else
        D = 0;

    // Kinda Stable
    count[1] = count[2] = 0;
                // std::ofstream log("TRY.txt");
    // for (int ply = 1; ply <= 2; ply++) {
    //     memset(visit,false,sizeof(visit));
    //     for (Point p : corners) {
    //         if (_board[p.x][p.y] == ply){
    //             count[ply] += countStableDisc(_board, p, ply);
    //             // log << count[ply]<< " ";
                
    //         }
    //     }
    // }
    if (count[curPlayer] > count[getNextPlayer(curPlayer)])
        SDC = (100.0 * count[curPlayer]) /
             (count[curPlayer] + count[getNextPlayer(curPlayer)]);
    else if (count[curPlayer] < count[getNextPlayer(curPlayer)])
        SDC = (100.0 * count[curPlayer]) /
             (count[curPlayer] + count[getNextPlayer(curPlayer)]);
    else
        SDC = 0;
    // Valid Moves Count
    count[curPlayer] = get_valid_spots(_board, curPlayer).size();
    count[getNextPlayer(curPlayer)] =
        get_valid_spots(_board, getNextPlayer(curPlayer)).size();
    if (count[curPlayer] > count[getNextPlayer(curPlayer)])
        MC = (100.0 * count[curPlayer]) /
            (count[curPlayer] + count[getNextPlayer(curPlayer)]);
    else if (count[curPlayer] < count[getNextPlayer(curPlayer)])
        MC = -(100.0 * count[curPlayer]) /
            (count[curPlayer] + count[getNextPlayer(curPlayer)]);
    else
        MC = 0;

    // Corner Instability
    count[1] = count[2] = 0;
    for (Point c : corners) {
        if (_board[c.x][c.y] == EMPTY) {
            for (int j = 1; j <= 7; j += 2) {
                Point p = c + directions[j];
                if (is_spot_on_board(p)) count[_board[p.x][p.y]]++;
            }
        }
    }

    CS = -20.25 * (count[curPlayer] - count[getNextPlayer(curPlayer)]);
    
    // Corners Captured
    count[1] = count[2] = 0;
    for (Point c : corners)
        count[_board[c.x][c.y]]++;


    C = 25 * (count[curPlayer] - count[getNextPlayer(curPlayer)]);

    double score = (10 * V) + (11 * D) + (80 * MC) + (375.78 * CS) + (805.131 * C) + (100.25 * SDC);
    return score;
}
Point StateValue()
{
    std::priority_queue<Node> pq;
    for (auto c : next_valid_spots) {
        ARR _board = board;
        put_disc(_board, c, player);
        pq.push(Node(c, Heuristic(_board, player)));
    }
    return pq.top().p;
}
double MiniMax(ARR _board, int depth, int curPlayer) {
    if (depth == 0) return Heuristic(_board, curPlayer);

    vector <Point> nextMove = get_valid_spots(_board, curPlayer);

    if (nextMove.size() == 0) return MiniMax(_board, depth - 1, getNextPlayer(curPlayer));
    //Maximizing
    if (curPlayer == player) {
        double val = -inf;
        for (Point p : nextMove) {
            ARR _state = _board;
            put_disc(_state, p, curPlayer);
            val = std::max(val, MiniMax(_state, depth - 1, getNextPlayer(curPlayer)));
        }
        return val;
    }
    // Minimizing
    else {
        double val = inf;
        for (Point p : nextMove) {
            ARR _state = _board;
            put_disc(_state, p , curPlayer);
            val = std::min(val, MiniMax(_state, depth - 1, getNextPlayer(curPlayer)));
        }
        return val;
    }
}
Point MiniMaxDecision(int depth, std::ofstream& fout) {
    double bestVal = -inf;

    for (auto p : next_valid_spots) {
        ARR _state = board;
        put_disc(_state, p , player);
        double tmp = MiniMax(_state, depth, getNextPlayer(player));

        if (tmp > bestVal) {
            bestMove = p;
            bestVal = tmp;
        }
        fout << bestMove.x << " " << bestMove.y << "\n";
    }
    return bestMove;
}
double AlphaBeta(ARR _board, int depth, int curPlayer, double a, double b) {
    if (depth == 0) return Heuristic(_board, curPlayer);

    vector <Point> nextMove = get_valid_spots(_board, curPlayer);
    if (nextMove.size() == 0) return AlphaBeta(_board, depth - 1, getNextPlayer(curPlayer), a, b);
    //Maximizing
    if (curPlayer == player) {
        double val = -inf;
        for (Point p : nextMove) {
            ARR _state = _board;
            put_disc(_state, p, curPlayer);
            val = std::max(val, AlphaBeta(_state, depth - 1, getNextPlayer(curPlayer), a, b));

            a = std::max(a, val);
            if (a >= b) break;
        }
        return val;
    }
    // Minimizing
    else {
        double val = inf;
        for (Point p : nextMove) {
            ARR _state = _board;
            put_disc(_state, p , curPlayer);
            val = std::min(val, AlphaBeta(_state, depth - 1, getNextPlayer(curPlayer), a, b));

            b = std::min(b, val);
            if (b <= a) break;
        }
        return val;
    }
}
Point AlphaBetaDecision(int depth, std::ofstream& fout) {
     double bestVal = -inf;

            // std::ofstream log("TRY.txt");
    for (auto p : next_valid_spots) {
        ARR _state = board;
        put_disc(_state, p , player);
        double tmp = AlphaBeta(_state, depth, getNextPlayer(player), -inf, inf);

        if (tmp > bestVal) {
            bestMove = p;
            bestVal = tmp;
        }

            // log << bestMove.x << " " << bestMove.y << "\n";

        fout << bestMove.x << " " << bestMove.y << "\n";
    }
    return bestMove;
}
/*----------------------------------------------------------------
double AlphaTabl(ARR _board, int depth, int curPlayer, double a, double b) {
    // Transposition Table Lookup
    if (Retrieve(_board)) {

    }
    
    if (depth == 0) return Heuristic(_board, curPlayer);

    vector <Point> nextMove = get_valid_spots(_board, curPlayer);
    if (nextMove.size() == 0) return AlphaBeta(_board, depth - 1, getNextPlayer(curPlayer), a, b);
    //Maximizing
    if (curPlayer == player) {
        double val = -inf;
        for (Point p : nextMove) {
            ARR _state = _board;
            put_disc(_state, p, curPlayer);
            val = std::max(val, AlphaBeta(_state, depth - 1, getNextPlayer(curPlayer), a, b));

            if (val >= b) break;
            a = std::max(a, val);
        }
        return val;
    }
    // Minimizing
    else {
        double val = inf;
        for (Point p : nextMove) {
            ARR _state = _board;
            put_disc(_state, p , curPlayer);
            val = std::min(val, AlphaBeta(_state, depth - 1, getNextPlayer(curPlayer), a, b));

            if (val <= a) break;
            b = std::min(b, val);
        }
        return val;
    }
}
double MTDF(int f, int depth, std::ofstream& fout) {
    // double lowerBound = -inf, upperBound = inf;
    // double g = d;
    // hash = computeBoardHash();
    // for (auto p : next_valid_spots) {
    //     ARR _state = board;
    //     put_disc(_state, p , player);
    //     double tmp = AlphaTabl(_state, depth, getNextPlayer(player), -inf, inf);

    //     if (tmp > bestVal) {
    //         bestMove = p;
    //         bestVal = tmp;
    //     }
    //     fout << bestMove.x << " " << bestMove.y << "\n";
    // }
    // return bestMove;
}
Point MTDFDecision(int depth, std::ostream& fout) {
    // double guess = 0;
    // Point bestMove
    // for (int d = 1; d <= depth; d++) {
    //     guess = MTDF
    // }
}
*/
void read_board(std::ifstream& fin)
{
    fin >> player;
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            fin >> board[i][j];
        }
    }
}

void read_valid_spots(std::ifstream& fin)
{
    int n_valid_spots;
    fin >> n_valid_spots;
    int x, y;
    for (int i = 0; i < n_valid_spots; i++) {
        fin >> x >> y;
        next_valid_spots.push_back({x, y});
    }
}

void write_valid_spot(std::ofstream& fout)
{
    // O is first player, X is second player
    Point p;
    if (algo == purerandom) {
        srand(time(NULL));
        int index = (rand() % next_valid_spots.size());
        p = next_valid_spots[index];
    }
    else if (algo == statevalue) {
        p = StateValue();
    }
    else if (algo == minimax) {
        p = MiniMaxDecision(DEPTH, fout);
        //  higher depths might run out of spaces 
    }
    else if (algo == alphabeta) {
        p = AlphaBetaDecision(DEPTH, fout);
    }
    else if (algo == mtdf) {
        // p = MTDFDecision(DEPTH, fout);
    }
    fout << p.x << " " << p.y << std::endl;
    fout.flush();
}
int main(int, char** argv)
{
    std::ifstream fin(argv[1]);
    std::ofstream fout(argv[2]);
    read_board(fin);
    read_valid_spots(fin);
    write_valid_spot(fout);
    fin.close();
    fout.close();
    return 0;
}
