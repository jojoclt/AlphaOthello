#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <queue>
#include <vector>

#define inf 1e9
const int SIZE = 8;
using ARR = std::array<std::array<int, SIZE>, SIZE>;
using std::vector;

enum Algo { purerandom, statevalue, minimax, alphabeta };
Algo algo = statevalue;

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
    Node() : p(0,0), val(0) {}
    Node(Point p, double val) : p(p), val(val) {}
};

int player = 1;
ARR board;
vector<Point> next_valid_spots;

enum SPOT_STATE { EMPTY = 0, BLACK = 1, WHITE = 2 };
const std::array<Point, 8> directions{{
        Point(-1, -1), Point(-1, 0), Point(-1, 1),
        Point(0, -1), /*{0, 0}, */Point(0, 1),
        Point(1, -1), Point(1, 0), Point(1, 1)
    }};

int get_next_player(int player) { return 3 - player; }
bool is_spot_on_board(Point p){ return 0 <= p.x && p.x < SIZE && 0 <= p.y && p.y < SIZE; }
int get_disc(ARR _board, Point p) { return _board[p.x][p.y]; }
void set_disc(ARR &_board, Point p, int disc) { _board[p.x][p.y] = disc; }
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
        if (!is_disc_at(_board, p, get_next_player(curPlayer)))
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
void flip_discs(ARR &_board, Point center, int curPlayer)
{
    for (Point dir : directions) {
        // Move along the direction while testing.
        Point p = center + dir;
        if (!is_disc_at(_board, p, get_next_player(curPlayer))) continue;
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
// Put disc on board and call flip
void put_disc(ARR& _board, Point p, int curPlayer)
{
    set_disc(_board, p, curPlayer);
    flip_discs(_board, p, curPlayer);
}
void print(ARR _board) {
    int i = 0;
   for (auto e : _board) {  
       for (auto c : e) {
           std::cout << c << " ";
       }
       printf("\n");
    }
    printf("=============================\n");
}
double Heuristic(ARR _board)
{
    int count[3] = {};
    double V = 0, D = 0, C = 0, S = 0, M = 0;
    std::array<std::array<int, SIZE>, SIZE> w;
    w[0] = {100, -10, 11, 6, 6, 11, -10, 100};
    w[1] = {-10, -20, -8, 2, 2, -8, -20, -10};
    w[2] = {10, 1, 5, 4, 4, 5, 1, 10};
    w[3] = {6, 2, 4, 2, 2, 4, 2, 6};
    w[4] = {6, 2, 4, 2, 2, 4, 2, 6};
    w[5] = {10, 1, 5, 4, 4, 5, 1, 10};
    w[6] = {-10, -20, -8, 2, 2, -8, -20, -10};
    w[7] = {100, -10, 11, 6, 6, 11, -10, 100};

    // Position Values and Pieces Count
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            if (_board[i][j] == player)
                V += w[i][j];
            else if (_board[i][j] == get_next_player(player))
                V -= w[i][j];
            count[_board[i][j]]++;
        }
    }
    if (count[player] > count[get_next_player(player)])
        D = (100.0 * count[player]) /
            (count[player] + count[get_next_player(player)]);
    else if (count[player] < count[get_next_player(player)])
        D = -(100.0 * count[player]) /
            (count[player] + count[get_next_player(player)]);
    else
        D = 0;

    // Valid Moves Count
    count[player] = get_valid_spots(_board, player).size();
    count[get_next_player(player)] =
        get_valid_spots(_board, get_next_player(player)).size();
    if (count[player] > count[get_next_player(player)])
        M = (100.0 * count[player]) /
            (count[player] + count[get_next_player(player)]);
    else if (count[player] < count[get_next_player(player)])
        M = -(100.0 * count[player]) /
            (count[player] + count[get_next_player(player)]);
    else
        M = 0;

    // Stability
    count[1] = count[2] = 0;
    if (_board[0][0] == EMPTY) {
        count[_board[0][1]]++;
        count[_board[1][1]]++;
        count[_board[1][0]]++;
    }
    if (_board[0][7] == EMPTY) {
        count[_board[0][6]]++;
        count[_board[1][6]]++;
        count[_board[1][7]]++;
    }
    if (_board[7][0] == EMPTY) {
        count[_board[7][1]]++;
        count[_board[6][1]]++;
        count[_board[6][0]]++;
    }
    if (_board[7][7] == EMPTY) {
        count[_board[6][7]]++;
        count[_board[6][6]]++;
        count[_board[7][6]]++;
    }
    S = -14.42 * (count[player] - count[get_next_player(player)]);
    // Corners Captured
    count[1] = count[2] = 0;
    count[_board[0][0]]++;
    count[_board[0][7]]++;
    count[_board[7][0]]++;
    count[_board[7][7]]++;

    C = count[player] - count[get_next_player(player)];

    double score = (10 * V) + (9 * D) + (88.98 * M) + (354.11 * S) + (1041.44 * C);
    return score;
}

double EvalMiniMax(ARR _state, int depth, int curPlayer) {
    if (depth == 0) return Heuristic(_state);

    vector<Point> nextMoves = get_valid_spots(_state, curPlayer);
    if (nextMoves.size() == 0)
        return EvalMiniMax(_state, depth - 1, get_next_player(curPlayer));

    double bestVal = (curPlayer == player) ? -inf : inf;
    for (auto c : nextMoves) {
        ARR _board = _state;
        put_disc(_board, c, curPlayer);
        double val = EvalMiniMax(_board, depth - 1, get_next_player(curPlayer));

        if (curPlayer == player)
            bestVal = std::max(bestVal, val);

        else 
            bestVal = std::min(bestVal, val);
    }
    return bestVal;
}
Point MiniMax(ARR board, int depth, int player) {
    Point bestMove;
    double bestVal = -inf;
    vector<Point> nextMoves = get_valid_spots(board, player);
   
    for (auto c : nextMoves) {
        ARR _board = board;
        put_disc(_board, c, player);
        double val = EvalMiniMax(_board, depth, player);

        if (val > bestVal) {
            bestVal = val;
            bestMove = c;
        }
        //  log << c.x << " " << c.y <<" " <<val << "\n";
    }
    // log << bestMove.x << " " << bestMove.y << "\n";
    std::cout<<"X";
    std::cout << bestMove.x << " " << bestMove.y <<"\n";
     return bestMove;
}

main() {
    ARR _board;
    _board = {0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 2, 0, 0,
              0, 0, 2, 2, 2, 1, 0, 0,
              0, 0, 0, 1, 1, 0, 0, 0,
              0, 0, 2, 1, 1, 0, 0, 0,
              0, 0, 0, 1, 0, 1, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0};
    // print(_board);
    MiniMax(_board, 4, 1);
}