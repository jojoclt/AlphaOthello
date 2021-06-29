#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <queue>
#include <vector>
#include<cassert>
#include<cstring>

#define inf 1e9

const int DEPTH = 3;
const int SIZE = 8;
using ARR = std::array<std::array<int, SIZE>, SIZE>;
using std::vector;

enum Algo { purerandom, statevalue, minimax, alphabeta, alphaenhanced };
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
int get_next_player(int player) { return 3 - player; }
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
void flip_discs(ARR& _board, Point center, int curPlayer)
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
void put_disc(ARR& _board, Point p, int curPlayer)
{
    set_disc(_board, p, curPlayer);
    flip_discs(_board, p, curPlayer);
}
bool visit[SIZE][SIZE] = {};
int CountStableDisc(ARR _board, Point p, int curPlayer) {
    std::queue <Point> q;
    q.push(p);
    int cnt = 0;
    while (!q.empty()) {
        Point p = q.front(); q.pop();
        if (visit[p.x][p.y]) continue;
        for (int i = 1; i <= 7; i += 2) {
            Point t = p + directions[i];
            if (!is_spot_on_board(t)) continue;
            if (visit[t.x][t.y]) continue;
            visit[t.x][t.y] = true;
            if (_board[t.x][t.y] != curPlayer) continue;
            ++cnt;
            q.push(t);
        }
    }
    return cnt;
}
double Heuristic(ARR _board, int curPlayer)
{
    memset(visit,false,sizeof(visit));
    int count[3] = {};
    int FS[3] = {};
    double V = 0, D = 0, C = 0, S = 0, M = 0, SS = 0, SC = 0;
    ARR w;
    w[0] = {20, -3, 11, 8, 8, 11, -3, 20};
    w[1] = {-3, -7, -4, -1, -1, -4, -7, -3};
    w[2] = {11, -4, 2, 2, 2, 2, -4, 11};
    w[3] = {8, 1, 2, -3, -3, 2, 1, 8};
    w[4] = {8, 1, 2, -3, -3, 2, 1, 8};
    w[5] = {11, -4, 2, 2, 2, 2, -4, 11};
    w[6] = {-3, -7, -4, -1, -1, -4, -7, -3};
    w[7] = {20, -3, 11, 8, 8, 11, -3, 20};

    // Position Values and Pieces Count
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            if (_board[i][j] == curPlayer)
                V += w[i][j];
            else if (_board[i][j] == get_next_player(curPlayer))
                V -= w[i][j];
            count[_board[i][j]]++;
            if(_board[i][j] != EMPTY)   {
				for(int k = 0; k < 8; k++)  {
					Point p = Point(i,j) + directions[k];
                    if (is_spot_on_board(p) && _board[i][j] == EMPTY) {
                        FS[_board[i][j]]++;
                        break;
                    }
				}
			}
        }
    }

    if (count[curPlayer] > count[get_next_player(curPlayer)])
        D = (100.0 * count[curPlayer]) /
            (count[curPlayer] + count[get_next_player(curPlayer)]);
    else if (count[curPlayer] < count[get_next_player(curPlayer)])
        D = -(100.0 * count[curPlayer]) /
            (count[curPlayer] + count[get_next_player(curPlayer)]);
    else
        D = 0;

    // Stability
    if (FS[curPlayer] > FS[get_next_player(curPlayer)])
        SS = -(100.0 * FS[curPlayer]) /
             (FS[curPlayer] + FS[get_next_player(curPlayer)]);
    else if (FS[curPlayer] < FS[get_next_player(curPlayer)])
        SS = (100.0 * FS[curPlayer]) /
             (FS[curPlayer] + FS[get_next_player(curPlayer)]);
    else
        SS = 0;

  // Count all stable disc
    count[1] = count[2] = 0;
    for (int i = 1; i <= 2; i++) {
        memset(visit,false,sizeof(visit));
        for (Point p : corners) {
            if (_board[p.x][p.y] == i){
                count[i] += CountStableDisc(_board, p, i);
            }
        }
    }
    if (count[curPlayer] > count[get_next_player(curPlayer)])
        SC = (100.0 * count[curPlayer]) /
             (count[curPlayer] + count[get_next_player(curPlayer)]);
    else if (count[curPlayer] < count[get_next_player(curPlayer)])
        SC = (100.0 * count[curPlayer]) /
             (count[curPlayer] + count[get_next_player(curPlayer)]);
    else
        SC = 0;

    // Valid Moves Count
    count[curPlayer] = get_valid_spots(_board, curPlayer).size();
    count[get_next_player(curPlayer)] =
        get_valid_spots(_board, get_next_player(curPlayer)).size();
    if (count[curPlayer] > count[get_next_player(curPlayer)])
        M = (100.0 * count[curPlayer]) /
            (count[curPlayer] + count[get_next_player(curPlayer)]);
    else if (count[curPlayer] < count[get_next_player(curPlayer)])
        M = -(100.0 * count[curPlayer]) /
            (count[curPlayer] + count[get_next_player(curPlayer)]);
    else
        M = 0;

    // Corner Stability
    count[1] = count[2] = 0;
    for (Point c : corners) {
        if (_board[c.x][c.y] == EMPTY) {
            for (int j = 1; j <= 7; j += 2) {
                Point p = c + directions[j];
                if (is_spot_on_board(p)) count[_board[p.x][p.y]]++;
            }
        }
    }

    S = -12.25 * (count[curPlayer] - count[get_next_player(curPlayer)]);
    
    // Corners Captured
    count[1] = count[2] = 0;
    for (Point c : corners)
        count[_board[c.x][c.y]]++;


    C = 25 * (count[curPlayer] - count[get_next_player(curPlayer)]);

    double score = (10 * V) + (10 * D) + (78.922 * M) + (382.026 * S) + (801.724 * C) + (74.396 * SS) + (301.25 * SC);
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

    if (nextMove.size() == 0) return MiniMax(_board, depth - 1, get_next_player(curPlayer));
    //Maximizing
    if (curPlayer == player) {
        double val = -inf;
        for (Point p : nextMove) {
            ARR _state = _board;
            put_disc(_state, p, curPlayer);
            val = std::max(val, MiniMax(_state, depth - 1, get_next_player(curPlayer)));
        }
        return val;
    }
    // Minimizing
    else {
        double val = inf;
        for (Point p : nextMove) {
            ARR _state = _board;
            put_disc(_state, p , curPlayer);
            val = std::min(val, MiniMax(_state, depth - 1, get_next_player(curPlayer)));
        }
        return val;
    }
}
Point MiniMaxDecision(int depth, std::ofstream& fout) {
    double bestVal = -inf;

    for (auto p : next_valid_spots) {
        ARR _state = board;
        put_disc(_state, p , player);
        double tmp = MiniMax(_state, depth, get_next_player(player));

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
    if (nextMove.size() == 0) return AlphaBeta(_board, depth - 1, get_next_player(curPlayer), a, b);
    //Maximizing
    if (curPlayer == player) {
        double val = -inf;
        for (Point p : nextMove) {
            ARR _state = _board;
            put_disc(_state, p, curPlayer);
            val = std::max(val, AlphaBeta(_state, depth - 1, get_next_player(curPlayer), a, b));

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
            val = std::min(val, AlphaBeta(_state, depth - 1, get_next_player(curPlayer), a, b));

            if (val <= a) break;
            b = std::min(b, val);
        }
        return val;
    }
}
Point AlphaBetaDecision(int depth, std::ofstream& fout) {
     double bestVal = -inf;

    for (auto p : next_valid_spots) {
        ARR _state = board;
        put_disc(_state, p , player);
        double tmp = AlphaBeta(_state, depth, get_next_player(player), -inf, inf);

        if (tmp > bestVal) {
            bestMove = p;
            bestVal = tmp;
        }
        fout << bestMove.x << " " << bestMove.y << "\n";
    }
    return bestMove;
}
double AlphaTabl(ARR _board, int depth, int curPlayer, double a, double b) {
    if (depth == 0) return Heuristic(_board, curPlayer);

    vector <Point> nextMove = get_valid_spots(_board, curPlayer);
    if (nextMove.size() == 0) return AlphaBeta(_board, depth - 1, get_next_player(curPlayer), a, b);
    //Maximizing
    if (curPlayer == player) {
        double val = -inf;
        for (Point p : nextMove) {
            ARR _state = _board;
            put_disc(_state, p, curPlayer);
            val = std::max(val, AlphaBeta(_state, depth - 1, get_next_player(curPlayer), a, b));

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
            val = std::min(val, AlphaBeta(_state, depth - 1, get_next_player(curPlayer), a, b));

            if (val <= a) break;
            b = std::min(b, val);
        }
        return val;
    }
}
Point AlphaTablDecision(int depth, std::ofstream& fout) {
     double bestVal = -inf;

    for (auto p : next_valid_spots) {
        ARR _state = board;
        put_disc(_state, p , player);
        double tmp = AlphaTabl(_state, depth, get_next_player(player), -inf, inf);

        if (tmp > bestVal) {
            bestMove = p;
            bestVal = tmp;
        }
        fout << bestMove.x << " " << bestMove.y << "\n";
    }
    return bestMove;
}

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
    int n_valid_spots = next_valid_spots.size();
    if (n_valid_spots == 0) return;

    Point p;
    if (algo == purerandom) {
        srand(time(NULL));
        int index = (rand() % n_valid_spots);
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
    else if (algo == alphaenhanced) {
        p = AlphaTablDecision(DEPTH, fout);
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
