#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <queue>
#include <vector>
#include<cassert>

#define inf 1e9
const int SIZE = 8;
using ARR = std::array<std::array<int, SIZE>, SIZE>;
using std::vector;

enum Algo { purerandom, statevalue, minimax, alphabeta };
Algo algo = minimax;

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
// Put disc on board and call flip
void put_disc(ARR& _board, Point p, int curPlayer)
{
    set_disc(_board, p, curPlayer);
    flip_discs(_board, p, curPlayer);
}
double Heuristic(ARR _board, int curPlayer)
{
    int count[3] = {};
    int FS[3] = {};
    double V = 0, D = 0, C = 0, S = 0, M = 0, SS = 0;
    std::array<std::array<int, SIZE>, SIZE> w;
    w[0] = {20, -3, 11, 8, 8, 11, -3, 20};
    w[1] = {-3, -7, -4, 1, 1, -4, -7, -3};
    w[2] = {11, -4, 2, 2, 2, 2, -4, 11};
    w[3] = {8, 1, 2, -3, -3, 2, 1, 8};
    w[4] = {8, 1, 2, -3, -3, 2, 1, 8};
    w[5] = {11, -4, 2, 2, 2, 2, -4, 11};
    w[6] = {-3, -7, -4, 1, 1, -4, -7, -3};
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
    S = -12.25 * (count[curPlayer] - count[get_next_player(curPlayer)]);
    // Corners Captured
    count[1] = count[2] = 0;
    count[_board[0][0]]++;
    count[_board[0][7]]++;
    count[_board[7][0]]++;
    count[_board[7][7]]++;

    C = 25 * (count[curPlayer] - count[get_next_player(curPlayer)]);

    double score = (10 * V) + (10 * D) + (78.922 * M) + (382.026 * S) + (801.724 * C) + (74.396 * SS);
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
    if (nextMove.size() == 0 && get_valid_spots(_board, get_next_player(curPlayer)).size() == 0) return 0;
    if (nextMove.size() == 0) return MiniMax(_board, depth - 1, get_next_player(curPlayer));
    //Maximizing
    if (curPlayer == player) {
        double val = -inf;
        for (Point p : nextMove) {
            ARR _state = _board;
            put_disc(_state, p, curPlayer);
            double eval = MiniMax(_state, depth - 1, get_next_player(curPlayer));
            if (eval > val) {
                val = eval;
            }
        }
        return val;
    }
    // Minimizing
    double val = inf;
    for (Point p : nextMove) {
        ARR _state = _board;
        put_disc(_state, p , curPlayer);
        double eval = MiniMax(_state, depth - 1, get_next_player(curPlayer));
        if (eval < val) {
            val = eval;
        }
    }
    return val;
}
Point MiniMaxDecision(int depth) {
    double bestVal = -inf;

    for (auto p : next_valid_spots) {
        ARR _state = board;
        put_disc(_state, p , player);
        double tmp = MiniMax(_state, depth, get_next_player(player));

        if (tmp > bestVal) {
            bestMove = p;
            bestVal = tmp;
        }
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
        p = MiniMaxDecision(3);
        //  higher depths might run out of spaces 
        //  std::ofstream log("TRY.txt");
        //  log << p.x << " " << p.y << "\n";
    }
    // else if (algo == alphabeta) {
    //     p = AlphaBeta(n_valid_spots);
    // }

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
