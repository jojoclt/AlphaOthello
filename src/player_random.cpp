#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <queue>
#include <vector>

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

struct node {
    Point p;
    double val;
    bool operator<(node t) const { return t.val > val; }
    node(Point p, double val) : p(p), val(val) {}
};

enum Algo { purerandom, statevalue, minimax, alphabeta };
Algo algo = statevalue;

int player;
const int SIZE = 8;
std::array<std::array<int, SIZE>, SIZE> board, _board;
std::vector<Point> next_valid_spots;

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
float Heuristic()
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
Point StateValue(int n)
{
    std::priority_queue<node> pq;

    for (int i = 0; i < n; i++) {
        _board = board;
        put_disc(next_valid_spots[i]);
        pq.push(node(next_valid_spots[i], Heuristic()));
    }
    return pq.top().p;
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

    Point p;
    if (algo == purerandom) {
        srand(time(NULL));
        int index = (rand() % n_valid_spots);
        p = next_valid_spots[index];
    }
    else if (algo == statevalue) {
        p = StateValue(n_valid_spots);
    }
    // else if (algo == minimax) {
    //     p = MiniMax(n_valid_spots);
    // }
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
