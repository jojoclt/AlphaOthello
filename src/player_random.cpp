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



int player;
ARR board, _board;
vector<Point> next_valid_spots;

enum SPOT_STATE { EMPTY = 0, BLACK = 1, WHITE = 2 };
const std::array<Point, 8> directions{{
        Point(-1, -1), Point(-1, 0), Point(-1, 1),
        Point(0, -1), /*{0, 0}, */Point(0, 1),
        Point(1, -1), Point(1, 0), Point(1, 1)
    }};
int get_next_player(int player) { return 3 - player; }
bool is_spot_on_board(Point p){ return 0 <= p.x && p.x < SIZE && 0 <= p.y && p.y < SIZE; }
int get_disc(Point p) { return _board[p.x][p.y]; }
void set_disc(Point p, int disc) { _board[p.x][p.y] = disc; }
bool is_disc_at(Point p, int disc)
{
    if (!is_spot_on_board(p)) return false;
    if (get_disc(p) != disc) return false;
    return true;
}
bool is_spot_valid(Point center, int curPlayer) {
    if (get_disc(center) != EMPTY)
        return false;
    for (Point dir: directions) {
        // Move along the direction while testing.
        Point p = center + dir;
        if (!is_disc_at(p, get_next_player(curPlayer)))
            continue;
        p = p + dir;
        while (is_spot_on_board(p) && get_disc(p) != EMPTY) {
            if (is_disc_at(p, curPlayer))
                return true;
            p = p + dir;
        }
    }
    return false;
}
int get_valid_spots_count(int curPlayer) {
    int cnt = 0;
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            Point p = Point(i, j);
            if (_board[i][j] != EMPTY)
                continue;
            if (is_spot_valid(p, curPlayer))
                cnt++;
        }
    }
    return cnt;
}
vector<Point> get_valid_spots(ARR _board, int curPlayer) {
    std::vector<Point> v;
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            Point p = Point(i, j);
            if (_board[i][j] != EMPTY)
                continue;
            if (is_spot_valid(p, curPlayer))
                v.push_back(p);
        }
    }
    return v;
}
void flip_discs(Point center)
{
    for (Point dir : directions) {
        // Move along the direction while testing.
        Point p = center + dir;
        if (!is_disc_at(p, get_next_player(player))) continue;
        std::vector<Point> discs({p});
        p = p + dir;
        while (is_spot_on_board(p) && get_disc(p) != EMPTY) {
            if (is_disc_at(p, player)) {
                for (Point s : discs) {
                    set_disc(s, player);
                }
                break;
            }
            discs.push_back(p);
            p = p + dir;
        }
    }
}
void put_disc(Point p)
{
    set_disc(p, player);
    flip_discs(p);
}
double Heuristic(ARR _board)
{
    int count[3] = {};
    double V = 0, D = 0, C = 0, M = 0;
    std::array<std::array<int, SIZE>, SIZE> w;
    w[0] = {100, -10, 11, 6, 6, 11, -10, 100};
    w[1] = {-10, -20, 1, 2, 2, 1, -20, -10};
    w[2] = {10, 1, 5, 4, 4, 5, 1, 10};
    w[3] = {6, 2, 4, 2, 2, 4, 2, 6};
    w[4] = {6, 2, 4, 2, 2, 4, 2, 6};
    w[5] = {10, 1, 5, 4, 4, 5, 1, 10};
    w[6] = {-10, -20, 1, 2, 2, 1, -20, -10};
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
    if(count[player] > count[get_next_player(player)])
		D = (100.0 * count[player])/(count[player] + count[get_next_player(player)]);
	else if(count[player] < count[get_next_player(player)])
		D = -(100.0 * count[player])/(count[player] + count[get_next_player(player)]);
	else D = 0;

    // Valid Moves Count
    count[player] = get_valid_spots_count(player);
    count[get_next_player(player)] = get_valid_spots_count(get_next_player(player));
    if (count[player] > count[get_next_player(player)])
        M = (100.0 * count[player])/(count[player] + count[get_next_player(player)]);
    else if (count[player] < count[get_next_player(player)])
        M = -(100.0 * count[player])/(count[player] + count[get_next_player(player)]);
    else M = 0;

    // Corners Captured
    count[1] = count[2] = 0;
    count[_board[0][0]]++;
    count[_board[0][7]]++;
    count[_board[7][0]]++;
    count[_board[7][7]]++;

    C = count[player] - count[get_next_player(player)];

    double score = (10 * V) + (10 * D) + (77.98 * M) + (752.44 * C);
    return score;
}
Point StateValue()
{
    std::priority_queue<Node> pq;

    for (int i = 0; i < next_valid_spots.size(); i++) {
        _board = board;
        put_disc(next_valid_spots[i]);
        pq.push(Node(next_valid_spots[i], Heuristic(_board)));
    }
    return pq.top().p;
}
double EvalMiniMax(ARR _state, int depth, int p) { 
    if (p == 0) {
        return Heuristic(_state);
    }
    vector<Point> nextMoves = get_valid_spots(_state, p);
    double bestVal = (p == player) ? -inf : inf;
    for (auto c : nextMoves) {
         _board = _state;
        put_disc(c);
        double val = EvalMiniMax(_board, depth - 1, get_next_player(p));

        if (p == player) 
            bestVal = std::max(bestVal, val);
        
        else 
            bestVal = std::min(bestVal, val);
    }

    return bestVal;
    
}
Point MiniMax(ARR _state, int depth, int p) {
    
    double bestVal = -inf;
    Point bestMove = next_valid_spots[0];
    for (auto c : next_valid_spots) {
        _board = _state;
        put_disc(c);

        double val = EvalMiniMax(_board, depth, p);

        if (val > bestVal) {
            bestVal = val;
            bestMove = c;
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
        p = MiniMax(board, 4, player);
    }
    // else if (algo == alphabeta) {
    //     p = AlphaBeta(n_valid_spots);
    // }
    // Remember to flush the output to ensure the last action is written to
    // file.
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
